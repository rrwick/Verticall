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
from .log import log, section_header, explanation
from .paint import paint_alignments, paint_assemblies


def pairwise(args):
    welcome_message(args)
    assemblies = find_assemblies(args.in_dir)
    build_indices(args, assemblies)
    with open(args.out_file, 'wt') as table_file:
        table_file.write(get_table_header())
        process_all_pairs(args, assemblies, table_file)
    finished_message()


def welcome_message(args):
    section_header('Starting Verticall pairwise')
    explanation('Vertical pairwise performs a pairwise analysis of all assemblies in the given '
                'directory, outputting the results to a tab-delimited table.')
    # TODO: make a different welcome message if --reference was used


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
        extension_len = len(extension) + 1  # plus one for the dot
        fasta_assemblies = sorted(f for f in in_dir.glob('*.' + extension) if f.is_file())
        for a in fasta_assemblies:
            sample_name = a.name[:-extension_len]
            if sample_name in all_assemblies:
                sys.exit(f'\nError: duplicate sample name {sample_name}')
            all_assemblies[sample_name] = a

    if extensions is None:
        extensions = ['fasta', 'fasta.gz', 'fna', 'fna.gz', 'fa', 'fa.gz']

    assemblies = {}
    for e in extensions:
        find_assemblies_with_extension(e, assemblies)
    assemblies = sorted(assemblies.items())

    log(f'Found {len(assemblies):,} samples in {in_dir.resolve()}\n')
    return assemblies


def process_all_pairs(args, assemblies, table_file):
    section_header('Processing pairwise combinations')
    explanation('For each assembly pair, Verticall pairwise aligns the assemblies, counts '
                'differences in a sliding window, builds a distribution and categories regions '
                'of the alignments as either vertical or horizontal. This allows for the '
                'calculation of a vertical-only genomic distance.')
    arg_list = []
    for name_a, filename_a in assemblies:
        for name_b, filename_b in assemblies:
            if name_a != name_b:
                arg_list.append((args, name_a, name_b, filename_a, filename_b))

    # If only using a single thread, do the alignment in a simple loop (easier for debugging).
    if args.threads == 1:
        for a in arg_list:
            log_text, name_a, name_b, distances, table_line = process_one_pair(a)
            log('\n'.join(prepare_log_text(log_text, args.verbose)), end='\n\n')
            table_file.write(table_line)
            table_file.flush()

    # If using multiple threads, use a process pool to work in parallel.
    else:
        with Pool(processes=args.threads) as pool:
            for log_text, name_a, name_b, distances, table_line in pool.imap(process_one_pair,
                                                                             arg_list):
                log('\n'.join(prepare_log_text(log_text, args.verbose)), end='\n\n')
                table_file.write(table_line)
                table_file.flush()


def prepare_log_text(log_text, verbose):
    prepared = []
    for line in log_text:
        if line.startswith('V'):
            if verbose:
                prepared.append(line[1:])
        else:
            prepared.append(line)
    return prepared


def process_one_pair(all_args, view=False):
    """
    This is the master function for each pairwise comparison. It gets called once for each assembly
    pair and carries out all analysis on that pair.

    Since this function can be called in parallel, it doesn't log text directly, but instead
    collects the text to be logged in a list which is returned.
    """
    args, name_a, name_b, filename_a, filename_b = all_args
    all_log_text = [f'{name_a} vs {name_b}:']

    alignments, n50_alignment_length, aligned_frac, log_text = \
        align_sample_pair(args, filename_a, name_b)
    all_log_text += log_text

    masses, window_size, window_count, mean_distance, median_distance, log_text = \
        get_distribution(args, alignments)
    all_log_text += log_text

    smoothed_masses = smooth_distribution(masses, args.smoothing_factor)
    mass_peaks, peak_distance, thresholds, log_text = \
        get_peak_distance(smoothed_masses, window_size)
    all_log_text += log_text

    vertical_masses, horizontal_masses, mean_vert_distance, median_vert_distance, log_text = \
        paint_alignments(alignments, thresholds, window_size)
    all_log_text += log_text

    painted_a, painted_b, log_text = \
        paint_assemblies(name_a, name_b, filename_a, filename_b, alignments)
    all_log_text += log_text

    # If being called by the view subcommand, we return the results instead of making a table line.
    if view:
        return alignments, window_size, masses, smoothed_masses, thresholds, vertical_masses, \
               horizontal_masses, painted_a, all_log_text

    distances = {'aligned_fraction': aligned_frac,
                 'mean': mean_distance,
                 'median': median_distance,
                 'peak': peak_distance,
                 'mean_vertical': mean_vert_distance,
                 'median_vertical': median_vert_distance}

    table_line = get_table_line(name_a, name_b, len(alignments), n50_alignment_length,
                                aligned_frac, window_size, window_count, mean_distance,
                                median_distance, mass_peaks, peak_distance, vertical_masses,
                                horizontal_masses, mean_vert_distance, median_vert_distance,
                                painted_a, painted_b)

    return all_log_text, name_a, name_b, distances, table_line


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
                   peak_distance, vertical_masses, horizontal_masses, mean_vert_distance,
                   median_vert_distance, painted_a, painted_b):
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
