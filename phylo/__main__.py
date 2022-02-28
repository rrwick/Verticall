#!/usr/bin/env python3
"""
Copyright 2022 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/XXXXXXXXX

This file is part of XXXXXXXXX. XXXXXXXXX is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. XXXXXXXXX is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with XXXXXXXXX.
If not, see <https://www.gnu.org/licenses/>.
"""

import argparse
from multiprocessing import Pool
import pathlib
import sys

from .align import build_indices
from .help_formatter import MyParser, MyHelpFormatter
from .misc import check_python_version, get_ascii_art, get_default_thread_count
from .log import bold, log, section_header, explanation
from .version import __version__


def main():
    check_python_version()
    args = parse_args(sys.argv[1:])
    check_args(args)
    welcome_message()
    assemblies = find_assemblies(args.in_dir)
    build_indices(args.in_dir, assemblies)
    process_all_pairs(args, assemblies)
    finished_message()


def parse_args(args):
    description = 'R|' + get_ascii_art() + '\n' + \
                  bold('XXXXXXXXX: a tool for finding recombination and generating '
                       'recombination-free trees')
    parser = MyParser(description=description, formatter_class=MyHelpFormatter, add_help=False)

    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument('-i', '--in_dir', type=pathlib.Path, required=True,
                               help='Directory containing assemblies in FASTA format')
    required_args.add_argument('-o', '--out_dir', type=pathlib.Path, required=True,
                               help='Output directory where results will be saved')

    align_args = parser.add_argument_group('Alignment settings')
    align_args.add_argument('--allowed_overlap', type=int, default=100,
                            help='Allow this much overlap between alignments')
    align_args.add_argument('--target_window_count', type=int, default=50000,
                            help='Aim to have at least this many comparison windows between '
                                 'assemblies')
    align_args.add_argument('--ignore_indels', action='store_true',
                            help='Only use mismatches to determine distance (default: use '
                                 'both mismatches and gap-compressed indels)')
    align_args.add_argument('--minimap2_options', type=str, default='TODO',
                            help='Minimap2 options for assembly-to-assembly alignment')

    distance_args = parser.add_argument_group('Distance settings')
    distance_args.add_argument('--method', type=str,
                               choices=['mean', 'median', 'mode', 'peak'], default='peak',
                               help='Method for converting distributions into a single distance')
    distance_args.add_argument('--correction', type=str, default='jukescantor',
                               help='Distance correction technique(s) from "none", "jukescantor" '
                                    'and "alignedfrac"')
    distance_args.add_argument('--asymmetrical', action='store_true',
                               help='Do not average pairs to make a symmetrical matrix')

    view_args = parser.add_argument_group('View settings')
    view_args.add_argument('--view', type=str,
                           help='Two assemblies (comma-separated) to view in plots')
    view_args.add_argument('--sqrt_x', action='store_true',
                           help='Use a square-root transform on the x-axis')
    view_args.add_argument('--sqrt_y', action='store_true',
                           help='Use a square-root transform on the y-axis')

    setting_args = parser.add_argument_group('General settings')
    setting_args.add_argument('-t', '--threads', type=int, default=get_default_thread_count(),
                              help='CPU threads for parallel processing')

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version', version='XXXXXXXXX v' + __version__,
                           help="Show program's version number and exit")

    # If no arguments were used, print the help.
    if len(args) == 0:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    return parser.parse_args(args)


def check_args(args):
    # Check --correction
    args.correction = set(args.correction.split(','))
    if 'none' in args.correction and len(args.correction) > 1:
        sys.exit('Error: --correction cannot contain "none" and additional values')
    valid_removed = set(args.correction)
    for c in ["none", "jukescantor", "alignedfrac"]:
        valid_removed.discard(c)
    if len(valid_removed) > 0:
        sys.exit('Error: only "none", "jukescantor" and "alignedfrac" can be used in '
                 '--correction')


def welcome_message():
    section_header('Starting XXXXXXXXX align')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')


def finished_message():
    section_header('Finished!')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')


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

    log(f'Found {len(assemblies):,} samples in {in_dir.resolve()}')
    log()
    return assemblies


def process_all_pairs(args, assemblies):
    section_header('Processing pairwise combinations')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')
    arg_list = []
    for sample_name_a, assembly_filename_a in assemblies:
        for sample_name_b, assembly_filename_b in assemblies:
            if sample_name_a != sample_name_b:
                arg_list.append((args, sample_name_a, sample_name_b))

    # If only using a single thread, do the alignment in a simple loop (easier for debugging).
    if args.threads == 1:
        for a in arg_list:
            log_text = process_one_pair(a)
            log('\n'.join(log_text), end='\n\n')

    # If using multiple threads, use a process pool to work in parallel.
    else:
        with Pool(processes=args.threads) as pool:
            for log_text in pool.imap(process_one_pair, arg_list):
                log('\n'.join(log_text), end='\n\n')


def process_one_pair(all_args):
    args, sample_name_a, sample_name_b = all_args
    log_text = [f'{sample_name_a} vs {sample_name_b}:']

    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO

    return log_text


if __name__ == '__main__':
    main()
