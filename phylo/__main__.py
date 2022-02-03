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
If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import pathlib
import sys

from .align import align
from .distance import distance
from .help_formatter import MyParser, MyHelpFormatter
from .misc import check_python_version, get_ascii_art
from .log import bold
from .version import __version__
from .view import view


def main():
    check_python_version()
    args = parse_args(sys.argv[1:])

    if args.subparser_name == 'align':
        check_align_args(args)
        align(args)

    if args.subparser_name == 'view':
        check_view_args(args)
        view(args)

    if args.subparser_name == 'distance':
        check_distance_args(args)
        distance(args)


def parse_args(args):
    description = 'R|' + get_ascii_art() + '\n' + \
                  bold('XXXXXXXXX: a tool for generating recombination-free trees')
    parser = MyParser(description=description, formatter_class=MyHelpFormatter, add_help=False)

    subparsers = parser.add_subparsers(title='Commands', dest='subparser_name')
    align_subparser(subparsers)
    view_subparser(subparsers)
    distance_subparser(subparsers)

    longest_choice_name = max(len(c) for c in subparsers.choices)
    subparsers.help = 'R|'
    for choice, choice_parser in subparsers.choices.items():
        padding = ' ' * (longest_choice_name - len(choice))
        subparsers.help += choice + ': ' + padding
        d = choice_parser.description
        subparsers.help += d[0].lower() + d[1:]  # don't capitalise the first letter
        subparsers.help += '\n'

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version', version='XXXXXXXXX v' + __version__,
                           help="Show program's version number and exit")

    # If no arguments were used, print the base-level help which lists possible commands.
    if len(args) == 0:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    return parser.parse_args(args)


def align_subparser(subparsers):
    group = subparsers.add_parser('align', description='align pieces to assemblies',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('-i', '--in_dir', type=pathlib.Path, required=True,
                               help='Directory containing assemblies in FASTA format')
    required_args.add_argument('-o', '--out_file', type=pathlib.Path, required=True,
                               help='Output file where alignment results will be saved')

    setting_args = group.add_argument_group('Settings')
    setting_args.add_argument('-t', '--threads', type=int, default=8,
                              help='Threads for alignment')
    setting_args.add_argument('--window_size', type=int, default=5000,
                              help='Size of the to sliding window over alignments')
    setting_args.add_argument('--window_step', type=int, default=1000,
                              help='Step distance of the sliding window over alignments')
    setting_args.add_argument('--min_align_len', type=int, default=1000,
                              help='Discard alignments shorter than this length')
    setting_args.add_argument('--allowed_overlap', type=int, default=100,
                              help='Allow this much overlap between alignments')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='XXXXXXXXX v' + __version__,
                            help="Show program's version number and exit")


def view_subparser(subparsers):
    group = subparsers.add_parser('view', description='view pairwise distance distribution',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('alignment_results', type=pathlib.Path,
                               help='File containing the output of XXXXXXXXX align')
    required_args.add_argument('assembly_1', type=str,
                               help='Name of first assembly in pair')
    required_args.add_argument('assembly_2', type=str,
                               help='Name of first assembly in pair')

    setting_args = group.add_argument_group('Settings')
    setting_args.add_argument('--sqrt_y', action='store_true',
                              help='Use a square-root transform on the y-axis')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='XXXXXXXXX v' + __version__,
                            help="Show program's version number and exit")


def distance_subparser(subparsers):
    group = subparsers.add_parser('distance', description='create PHYLIP distance matrix from '
                                                          'alignment results',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('alignment_results', type=pathlib.Path,
                               help='File containing the output of XXXXXXXXX align')

    setting_args = group.add_argument_group('Settings')
    setting_args.add_argument('--method', type=str, choices=['mean', 'median', 'median_cont'],
                              default='median_cont',
                              help='Method for converting distributions into a single distance')
    setting_args.add_argument('--correction', type=str, choices=['none', 'jukescantor'],
                              default='jukescantor',
                              help='Distance correction technique')
    setting_args.add_argument('--asymmetrical', action='store_true',
                              help='Do not average distance pairs to make a symmetrical matrix')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='XXXXXXXXX v' + __version__,
                            help="Show program's version number and exit")


def check_align_args(args):
    pass


def check_view_args(args):
    if not args.alignment_results.is_file():
        sys.exit(f'\nError: {args.alignment_results} could not be found')


def check_distance_args(args):
    if not args.alignment_results.is_file():
        sys.exit(f'\nError: {args.alignment_results} could not be found')


if __name__ == '__main__':
    main()
