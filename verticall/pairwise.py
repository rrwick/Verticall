#!/usr/bin/env python3
"""
This module contains the code for the 'verticall pairwise' subcommand.

Copyright 2022 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Verticall

This file is part of Verticall. Verticall is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Verticall is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Verticall.
If not, see <https://www.gnu.org/licenses/>.
"""

from multiprocessing import Pool
import sys

from .alignment import build_indices, align_sample_pair
from .distance import get_distribution, smooth_distribution, get_peak_distance
from .log import log, section_header, explanation, warning
from .misc import split_list, iterate_fasta, contains_ambiguous_bases
from .paint import paint_alignments, paint_assemblies


def pairwise(args):
    welcome_message(args)
    assemblies = find_assemblies(args.in_dir)
    reference = find_reference(args.reference)
    check_assemblies(assemblies, reference)
    build_indices(args, assemblies)
    with open(args.out_file, 'wt') as table_file:
        if parse_part(args.part)[0] == 0:  # only include the header in the first part
            table_file.write(get_table_header())
        process_all_pairs(args, assemblies, reference, table_file)
    finished_message()


def welcome_message(args):
    section_header('Starting Verticall pairwise')
    if args.reference is None:
        explanation_text = 'Vertical pairwise performs a pairwise analysis of all assemblies in ' \
                           'the given directory, outputting the results to a tab-delimited table.'
    else:
        explanation_text = 'Vertical pairwise performs a pairwise analysis of each assembly in ' \
                           'the given directory to the specified reference genome, outputting ' \
                           'the results to a tab-delimited table.'

    part_num, part_total = parse_part(args.part)
    if part_total > 1:
        explanation_text += f' Because --part {args.part} was used, this command will only ' \
                            f'analyse a fraction of the pairwise combinations.'

    explanation(explanation_text)


def finished_message():
    section_header('Finished!')
    explanation('You can now use the resulting tab-delimited file to produce a distance matrix '
                '(using Vertical matrix), summarise an assembly\'s vertical-vs-horizontal regions '
                '(using Vertical summary) or mask horizontal regions from a SNV matrix (using '
                'Vertical mask).')


def find_assemblies(in_dir, extensions=None):
    """
    Returns assemblies in a (sample_name, filename) tuple.
    """
    def find_assemblies_with_extension(extension, all_assemblies):
        fasta_assemblies = sorted(f for f in in_dir.glob('*' + extension) if f.is_file())
        for a in fasta_assemblies:
            sample_name = a.name[:-len(extension)]
            if sample_name in all_assemblies:
                sys.exit(f'Error: duplicate sample name {sample_name}')
            all_assemblies[sample_name] = a

    if extensions is None:
        extensions = get_default_assembly_extensions()

    assemblies = {}
    for e in extensions:
        find_assemblies_with_extension(e, assemblies)
    assemblies = sorted(assemblies.items())

    log(f'Found {len(assemblies):,} samples in {in_dir.resolve()}\n')
    return assemblies


def find_reference(reference, extensions=None):
    if reference is None:
        return None
    if extensions is None:
        extensions = get_default_assembly_extensions()
    for extension in extensions:
        if str(reference).endswith(extension):
            sample_name = reference.name[:-len(extension)]
            return sample_name, reference
    sys.exit(f'Error: {reference} does not end in a FASTA file extension')


def get_default_assembly_extensions():
    return ['.fasta', '.fasta.gz', '.fna', '.fna.gz', '.fa', '.fa.gz']


def check_assemblies(assemblies, reference=None):
    """
    Checks to make sure the assemblies look good: no duplicate contig names, no ambiguous bases.
    """
    if reference is not None:
        assemblies = sorted(set(assemblies + [reference]))
    files_with_duplicate_contig_names, files_with_ambiguous_bases = [], []
    log(f'Checking assemblies: 0 / {len(assemblies)}', end='')
    for i, a in enumerate(assemblies):
        sample_name, filename = a
        duplicate_contig_names, ambiguous_bases = check_one_assembly(filename)
        if duplicate_contig_names:
            files_with_duplicate_contig_names.append(filename)
        if ambiguous_bases:
            files_with_ambiguous_bases.append(filename)
        log(f'\rChecking assemblies: {i + 1} / {len(assemblies)}', end='')
    log('\n')
    if not files_with_duplicate_contig_names and not files_with_ambiguous_bases:
        return
    error_message = []
    if files_with_duplicate_contig_names:
        file_list = ', '.join([f.name for f in files_with_duplicate_contig_names])
        error_message.append(f'The following assemblies have duplicate contig names: {file_list}\n')
    if files_with_ambiguous_bases:
        file_list = ', '.join([f.name for f in files_with_ambiguous_bases])
        error_message.append(f'The following assemblies have ambiguous bases: {file_list}\n')
    error_message.append('Error: You must run verticall repair to fix assembly problems before '
                         'these assemblies can be used.')
    sys.exit('\n'.join(error_message))


def check_one_assembly(filename):
    """
    Returns a tuple of two booleans:
    * True if there are duplicates contig names, False if the contig names are okay.
    * True if there are ambiguous bases, False if the sequences are okay.
    """
    contig_names = []
    ambiguous_bases = False
    for name, seq in iterate_fasta(filename):
        contig_names.append(name)
        if ambiguous_bases or contains_ambiguous_bases(seq):
            ambiguous_bases = True
    duplicate_contig_names = len(contig_names) > len(set(contig_names))
    return duplicate_contig_names, ambiguous_bases


def process_all_pairs(args, assemblies, reference, table_file):
    section_header('Processing pairwise combinations')
    explanation('For each assembly pair, Verticall pairwise aligns the assemblies, counts '
                'differences in a sliding window, builds a distribution and categorises regions '
                'of the alignments as either vertical or horizontal. This allows for the '
                'calculation of a vertical-only genomic distance.')
    arg_list = get_arg_list(args, assemblies, reference)
    empty_results, multi_results = False, False

    # If only using a single thread, do the alignment in a simple loop (easier for debugging).
    if args.threads == 1:
        for a in arg_list:
            log_text, table_lines = process_one_pair(a)
            if len(table_lines) == 0:
                empty_results = True
            if len(table_lines) > 1:
                multi_results = True
            log('\n'.join(prepare_log_text(log_text, args.verbose)), end='\n\n')
            for table_line in table_lines:
                table_file.write(table_line)
            table_file.flush()

    # If using multiple threads, use a process pool to work in parallel.
    else:
        with Pool(processes=args.threads) as pool:
            for log_text, table_lines in pool.imap(process_one_pair, arg_list):
                if len(table_lines) == 0:
                    empty_results = True
                if len(table_lines) > 1:
                    multi_results = True
                log('\n'.join(prepare_log_text(log_text, args.verbose)), end='\n\n')
                for table_line in table_lines:
                    table_file.write(table_line)
                table_file.flush()

    if empty_results:
        warning('one or more assembly pairs failed to align sufficiently to produce results')
    if multi_results:
        warning('one or more assembly pairs produced multiple results')


def get_arg_list(args, assemblies, reference):
    """
    This function produces a list of arguments for the process_one_pair function. If --part 1/1 was
    used (the default), this will include an entry for each pair of assemblies. If another value
    for --part was used, this will include a subset of the pairs.
    """
    arg_list = []
    if reference is None:
        for name_a, filename_a in assemblies:
            for name_b, filename_b in assemblies:
                if name_a != name_b:
                    arg_list.append((args, name_a, name_b, filename_a, filename_b))
    else:
        ref_name, ref_filename = reference
        for assembly_name, assembly_filename in assemblies:
            if ref_name != assembly_name:
                arg_list.append((args, ref_name, assembly_name, ref_filename, assembly_filename))

    part_num, part_total = parse_part(args.part)
    if part_total > 1:
        arg_list = split_list(arg_list, part_total)[part_num]

    return arg_list


def parse_part(part_str):
    """
    Returns the numerator and denominator from the --part argument. The numerator is returned as a
    zero-based index, e.g. '1/1' -> (0, 1).
    """
    try:
        numerator, denominator = part_str.split('/')
        numerator, denominator = int(numerator)-1, int(denominator)
    except ValueError:
        sys.exit('Error: --part must be a fraction, e.g. 1/1, 1/10, 3/10, etc.')
    if numerator < 0:
        sys.exit('Error: the numerator of --part must be a positive integer')
    if denominator < 1:
        sys.exit('Error: the denominator of --part must be a positive integer')
    if numerator >= denominator:
        sys.exit('Error: the numerator of --part must be less than or equal to the denominator')
    return numerator, denominator


def prepare_log_text(log_text, verbose):
    prepared = []
    for line in log_text:
        if line.startswith('V'):
            if verbose:
                prepared.append(line[1:])
        else:
            prepared.append(line)
    return prepared


def process_one_pair(all_args, view=False, view_num=1):
    """
    This is the master function for each pairwise comparison. It gets called once for each assembly
    pair and carries out all analysis on that pair.

    Since this function can be run in parallel, it doesn't log text directly, but instead collects
    the text to be logged in a list which is returned.

    This function returns differently if it was called by the pairwise or view subcommands:
    * pairwise: returns a list of log text and a list of tsv table lines
    * view: returns the analysis data (alignments, mass distribution, painted assemblies, etc) to
            be plotted
    """
    args, name_a, name_b, filename_a, filename_b = all_args  # unpack the arguments
    all_log_text = [f'{name_a} vs {name_b}:']  # text logged to console will be stored in this list

    # Step 1: align the two assemblies to each other.
    alignments, n50_alignment_length, aligned_frac, log_text = \
        align_sample_pair(args, filename_a, name_b)
    all_log_text += log_text

    # Step 2: produce a distance distribution from sliding windows across the alignments.
    masses, window_size, window_count, mean_distance, median_distance, log_text = \
        get_distribution(args, alignments)
    all_log_text += log_text

    # If there are no sliding windows, then the two assemblies didn't sufficiently align to do any
    # further analysis.
    if window_count == 0:
        return all_log_text, []

    # Step 3: smooth the distribution and find peaks with their corresponding thresholds. When
    #         there is a close call, this can return multiple results (a primary result and one or
    #         more secondary results).
    smoothed_masses = smooth_distribution(masses, args.smoothing_factor)
    mass_peaks, results, log_text = get_peak_distance(smoothed_masses, window_size, args.secondary)
    all_log_text += log_text
    check_view_num(view, view_num, len(results))

    # Step 4: paint alignments and assemblies using the distance thresholds.
    table_lines = []
    for i, result in enumerate(results):
        if view and i+1 != view_num:
            continue

        _, result_level, peak_distance, thresholds = result

        vertical_masses, horizontal_masses, mean_vert_distance, median_vert_distance, log_text = \
            paint_alignments(alignments, thresholds, window_size)
        if result_level == 'primary' or args.verbose or view:
            all_log_text += log_text

        painted_a, painted_b, log_text = \
            paint_assemblies(name_a, name_b, filename_a, filename_b, alignments)
        if result_level == 'primary' or args.verbose or view:
            all_log_text += log_text

        # If called by the view subcommand, we return the results instead of making a table line.
        if view:
            return alignments, window_size, masses, smoothed_masses, thresholds, vertical_masses, \
                   horizontal_masses, painted_a, all_log_text

        table_lines.append(get_table_line(name_a, name_b, len(alignments), n50_alignment_length,
                                          aligned_frac, window_size, window_count, mean_distance,
                                          median_distance, mass_peaks, result_level, peak_distance,
                                          vertical_masses, horizontal_masses, mean_vert_distance,
                                          median_vert_distance, painted_a, painted_b))

    return all_log_text, table_lines


def check_view_num(view, view_num, result_count):
    if not view:
        return
    if result_count == 0:
        sys.exit(f'Error: this pair has no results and cannot be viewed')
    if view_num > result_count:
        result_str = 'result' if result_count == 1 else 'results'
        sys.exit(f'Error: this pair only has {result_count} {result_str} but --result {view_num} '
                 f'was used')


def get_table_header():
    return ('assembly_a\t'
            'assembly_b\t'
            'alignment_count\t'
            'n50_alignment_length\t'
            'aligned_fraction\t'
            'window_size\t'
            'window_count\t'
            'mean_distance\t'
            'median_distance\t'
            'mass_peaks\t'
            'result_level\t'
            'peak_distance\t'
            'alignments_vertical_fraction\t'
            'alignments_horizontal_fraction\t'
            'mean_vertical_distance\t'
            'median_vertical_distance\t'
            'assembly_a_vertical_fraction\t'
            'assembly_a_horizontal_fraction\t'
            'assembly_a_unaligned_fraction\t'
            'assembly_b_vertical_fraction\t'
            'assembly_b_horizontal_fraction\t'
            'assembly_b_unaligned_fraction\t'
            'assembly_a_vertical_regions\t'
            'assembly_a_horizontal_regions\t'
            'assembly_a_unaligned_regions\t'
            'assembly_b_vertical_regions\t'
            'assembly_b_horizontal_regions\t'
            'assembly_b_unaligned_regions\n')


def get_table_line(name_a, name_b, alignment_count, n50_alignment_length, aligned_frac,
                   window_size, window_count, mean_distance, median_distance, mass_peaks,
                   result_level, peak_distance, vertical_masses, horizontal_masses,
                   mean_vert_distance, median_vert_distance, painted_a, painted_b):
    vertical_frac_a, horizontal_frac_a, unaligned_frac_a = painted_a.get_fractions()
    vertical_frac_b, horizontal_frac_b, unaligned_frac_b = painted_b.get_fractions()
    vertical_regions_a, horizontal_regions_a, unaligned_regions_a = painted_a.get_regions()
    vertical_regions_b, horizontal_regions_b, unaligned_regions_b = painted_b.get_regions()
    return (f'{name_a}\t'
            f'{name_b}\t'
            f'{alignment_count}\t'
            f'{n50_alignment_length}\t'
            f'{aligned_frac:.9f}\t'
            f'{window_size}\t'
            f'{window_count}\t'
            f'{mean_distance:.9f}\t'
            f'{median_distance:.9f}\t'
            f'{mass_peaks}\t'
            f'{result_level}\t'
            f'{peak_distance:.9f}\t'
            f'{100.0 * sum(vertical_masses):.2f}%\t'
            f'{100.0 * sum(horizontal_masses):.2f}%\t'
            f'{mean_vert_distance:.9f}\t'
            f'{median_vert_distance:.9f}\t'
            f'{100.0 * vertical_frac_a:.2f}%\t'
            f'{100.0 * horizontal_frac_a:.2f}%\t'
            f'{100.0 * unaligned_frac_a:.2f}%\t'
            f'{100.0 * vertical_frac_b:.2f}%\t'
            f'{100.0 * horizontal_frac_b:.2f}%\t'
            f'{100.0 * unaligned_frac_b:.2f}%\t'
            f'{vertical_regions_a}\t'
            f'{horizontal_regions_a}\t'
            f'{unaligned_regions_a}\t'
            f'{vertical_regions_b}\t'
            f'{horizontal_regions_b}\t'
            f'{unaligned_regions_b}\n')
