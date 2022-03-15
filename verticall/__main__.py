#!/usr/bin/env python3
"""
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

import argparse
from multiprocessing import Pool
import pathlib
import sys

from .alignment import build_indices, align_sample_pair
from .distance import get_distribution, smooth_distribution, get_peak_distance
from .help_formatter import MyParser, MyHelpFormatter
from .matrix import save_all_matrices
from .misc import check_python_version, get_ascii_art, get_default_thread_count
from .log import bold, log, section_header, explanation
from .paint import paint_alignments, paint_assemblies
from .version import __version__
from .view import show_plots


def main():
    check_python_version()
    args = parse_args(sys.argv[1:])
    check_args(args)
    welcome_message()
    assemblies = find_assemblies(args.in_dir)
    if not args.view:
        check_output_directory(args.out_dir)
        table_file = create_table_file(args.out_dir)
    else:
        table_file = None
    build_indices(args, assemblies)
    if args.view is not None:
        view_one_pair(args, assemblies)
    else:
        sample_names, distances = process_all_pairs(args, assemblies, table_file)
        table_file.close()
        save_all_matrices(args, sample_names, distances)
    finished_message()


def parse_args(args):
    description = 'R|' + get_ascii_art() + '\n' + \
                  bold('Verticall: a tool for finding distance by vertical inheritance')
    parser = MyParser(description=description, formatter_class=MyHelpFormatter, add_help=False)

    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument('-i', '--in_dir', type=pathlib.Path, required=True,
                               help='Directory containing assemblies in FASTA format')
    required_args.add_argument('-o', '--out_dir', type=pathlib.Path, required=True,
                               help='Output directory where results will be saved')

    align_args = parser.add_argument_group('Alignment settings')
    # TODO: explore different indexing options (e.g. -k and -w) to see how they affect the results.
    align_args.add_argument('--index_options', type=str, default='-k15 -w10',
                            help='Minimap2 options for assembly indexing')
    # TODO: explore different alignment options (e.g. the things set by -x asm20).
    align_args.add_argument('--align_options', type=str, default='-x asm20',
                            help='Minimap2 options for assembly-to-assembly alignment')
    align_args.add_argument('--allowed_overlap', type=int, default=100,
                            help='Allow this much overlap between alignments')
    align_args.add_argument('--window_count', type=int, default=50000,
                            help='Aim to have at least this many comparison windows between '
                                 'assemblies')
    align_args.add_argument('--ignore_indels', action='store_true',
                            help='Only use mismatches to determine distance (default: use '
                                 'both mismatches and gap-compressed indels)')

    distance_args = parser.add_argument_group('Distance settings')
    distance_args.add_argument('--smoothing_factor', type=float, default=0.8,
                               help='Degree to which the distance distribution is smoothed')
    distance_args.add_argument('--asymmetrical', action='store_true',
                               help='Do not average pairs to make symmetrical matrices (default: '
                                    'make matrices symmetrical)')

    view_args = parser.add_argument_group('View settings')
    view_args.add_argument('--view', type=str,
                           help='Two assemblies (comma-separated) to view in plots')
    view_args.add_argument('--sqrt_distance', action='store_true',
                           help='Use a square-root transform on the genomic distance axis '
                                '(default: no distance axis transform)')
    view_args.add_argument('--sqrt_mass', action='store_true',
                           help='Use a square-root transform on the probability mass axis '
                                '(default: no mass axis transform)')

    setting_args = parser.add_argument_group('General settings')
    setting_args.add_argument('--verbose', action='store_true',
                              help='Output more detail to stderr for debugging (default: only '
                                   'output essential information)')
    setting_args.add_argument('-t', '--threads', type=int, default=get_default_thread_count(),
                              help='CPU threads for parallel processing')

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version', version='Verticall v' + __version__,
                           help="Show program's version number and exit")

    # If no arguments were used, print the help.
    if len(args) == 0:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    return parser.parse_args(args)


def check_args(args):
    if args.view is not None:
        samples = args.view.split(',')
        if len(samples) != 2:
            sys.exit('Error: two assemblies (comma-delimited) must be supplied to --view')


def welcome_message():
    section_header('Starting Verticall')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')


def finished_message():
    section_header('Finished!')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')


def check_output_directory(directory: pathlib.Path):
    if directory.is_file():
        sys.exit(f'\nError: output directory ({directory}) already exists as a file')
    if directory.is_dir():
        if len(list(directory.iterdir())) > 0:
            log(f'Output directory ({directory}) already exists and is not empty - files may be '
                f'overwritten.')
        else:
            log(f'Output directory already exists: {directory}')
    else:
        log(f'Creating output directory: {directory}')
        directory.mkdir(parents=True)
    log()


def create_table_file(directory: pathlib.Path):
    f = open(directory / 'pairwise.tsv', 'wt')
    f.write(get_table_header())
    return f


def find_assemblies(in_dir):
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

    assemblies = {}
    find_assemblies_with_extension('fasta', assemblies)
    find_assemblies_with_extension('fasta.gz', assemblies)
    find_assemblies_with_extension('fna', assemblies)
    find_assemblies_with_extension('fna.gz', assemblies)
    assemblies = sorted(assemblies.items())

    log(f'Found {len(assemblies):,} samples in {in_dir.resolve()}\n')
    return assemblies


def process_all_pairs(args, assemblies, table_file):
    section_header('Processing pairwise combinations')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')
    arg_list, sample_names = [], set()
    for name_a, filename_a in assemblies:
        sample_names.add(name_a)
        for name_b, filename_b in assemblies:
            if name_a != name_b:
                arg_list.append((args, name_a, name_b, filename_a, filename_b))

    all_distances = {}

    # If only using a single thread, do the alignment in a simple loop (easier for debugging).
    if args.threads == 1:
        for a in arg_list:
            log_text, name_a, name_b, distances, table_line = process_one_pair(a)
            log('\n'.join(prepare_log_text(log_text, args.verbose)), end='\n\n')
            all_distances[(name_a, name_b)] = distances
            table_file.write(table_line)
            table_file.flush()

    # If using multiple threads, use a process pool to work in parallel.
    else:
        with Pool(processes=args.threads) as pool:
            for log_text, name_a, name_b, distances, table_line in pool.imap(process_one_pair,
                                                                             arg_list):
                log('\n'.join(prepare_log_text(log_text, args.verbose)), end='\n\n')
                all_distances[(name_a, name_b)] = distances
                table_file.write(table_line)
                table_file.flush()

    return sorted(sample_names), all_distances


def prepare_log_text(log_text, verbose):
    prepared = []
    for line in log_text:
        if line.startswith('V'):
            if verbose:
                prepared.append(line[1:])
        else:
            prepared.append(line)
    return prepared


def view_one_pair(args, assemblies):
    section_header('Viewing single pair')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')
    name_a, name_b = args.view.split(',')
    filename_a = [filename for name, filename in assemblies if name == name_a][0]
    filename_b = [filename for name, filename in assemblies if name == name_b][0]
    process_one_pair([args, name_a, name_b, filename_a, filename_b], view=True)


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

    # TODO: get mean distance using vertical regions
    # TODO: save painting info to file

    if view:
        log('\n'.join(prepare_log_text(all_log_text, True)), end='\n\n')
        show_plots(name_a, name_b, window_size, aligned_frac, masses, smoothed_masses, thresholds,
                   vertical_masses, horizontal_masses, painted_a, painted_b, args.sqrt_distance,
                   args.sqrt_mass)

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
            'alignments_vertical\t'
            'alignments_horizontal\t'
            'mean_vertical_distance\t'
            'median_vertical_distance\t'
            'assembly_a_vertical\t'
            'assembly_a_horizontal\t'
            'assembly_a_unaligned\t'
            'assembly_b_vertical\t'
            'assembly_b_horizontal\t'
            'assembly_b_unaligned\n')


def get_table_line(name_a, name_b, alignment_count, n50_alignment_length, aligned_frac,
                   window_size, window_count, mean_distance, median_distance, mass_peaks,
                   peak_distance, vertical_masses, horizontal_masses, mean_vert_distance,
                   median_vert_distance, painted_a, painted_b):
    vertical_a, horizontal_a, unaligned_a = painted_a.get_fractions()
    vertical_b, horizontal_b, unaligned_b = painted_b.get_fractions()
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
            f'{100.0 * vertical_a:.2f}%\t'
            f'{100.0 * horizontal_a:.2f}%\t'
            f'{100.0 * unaligned_a:.2f}%\t'
            f'{100.0 * vertical_b:.2f}%\t'
            f'{100.0 * horizontal_b:.2f}%\t'
            f'{100.0 * unaligned_b:.2f}%\n')


if __name__ == '__main__':
    main()
