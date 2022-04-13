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
import pathlib
import sys

from .help_formatter import MyParser, MyHelpFormatter
from .misc import check_python_version, get_ascii_art, get_default_thread_count
from .log import bold
from .matrix import matrix
from .pairwise import pairwise
from .repair import repair
from .version import __version__
from .view import view, check_hex_colour


def main():
    check_python_version()
    args = parse_args(sys.argv[1:])

    if args.subparser_name == 'pairwise':
        check_pairwise_args(args)
        pairwise(args)
    elif args.subparser_name == 'view':
        check_view_args(args)
        view(args)
    elif args.subparser_name == 'matrix':
        check_matrix_args(args)
        matrix(args)
    # elif args.subparser_name == 'mask':
    #     check_mask_args(args)
    #     filter(args)
    # elif args.subparser_name == 'summary':
    #     check_summary_args(args)
    #     summary(args)
    elif args.subparser_name == 'repair':
        check_repair_args(args)
        repair(args)


def parse_args(args):
    description = 'R|' + get_ascii_art() + '\n' + \
                  bold('Verticall: a tool for finding vertical inheritance in bacterial genomes')
    parser = MyParser(description=description, formatter_class=MyHelpFormatter, add_help=False)

    subparsers = parser.add_subparsers(title='Commands', dest='subparser_name')
    pairwise_subparser(subparsers)
    view_subparser(subparsers)
    matrix_subparser(subparsers)
    mask_subparser(subparsers)
    summary_subparser(subparsers)
    repair_subparser(subparsers)

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
    help_args.add_argument('--version', action='version', version='Verticall v' + __version__,
                           help="Show program's version number and exit")

    # If no arguments were used, print the base-level help which lists possible commands.
    if len(args) == 0:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    return parser.parse_args(args)


def pairwise_subparser(subparsers):
    group = subparsers.add_parser('pairwise', description='pairwise analysis of assemblies',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('-i', '--in_dir', type=pathlib.Path, required=True,
                               help='Directory containing assemblies in FASTA format')
    required_args.add_argument('-o', '--out_file', type=pathlib.Path, required=True,
                               help='Filename for tsv output')

    reference_args = group.add_argument_group('Reference-based analysis')
    reference_args.add_argument('-r', '--reference', type=pathlib.Path,
                                help='Reference assembly in FASTA format')

    pairwise_and_view_settings(group)

    performance_args = group.add_argument_group('Performance')
    performance_args.add_argument('-t', '--threads', type=int, default=get_default_thread_count(),
                                  help='CPU threads for parallel processing')
    performance_args.add_argument('--part', type=str, default='1/1',
                                  help='Fraction of the data to analyse (for parallelisation, '
                                       'default: DEFAULT)')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='Verticall v' + __version__,
                            help="Show program's version number and exit")


def view_subparser(subparsers):
    group = subparsers.add_parser('view', description='view plots for a single assembly pair',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('-i', '--in_dir', type=pathlib.Path, required=True,
                               help='Directory containing assemblies in FASTA format')
    required_args.add_argument('-n', '--names', type=str, required=True,
                               help='Two sample names (comma-delimited) to be viewed')

    pairwise_and_view_settings(group)

    view_args = group.add_argument_group('Plot settings')
    view_args.add_argument('--sqrt_distance', action='store_true',
                           help='Use a square-root transform on the genomic distance axis '
                                '(default: no distance axis transform)')
    view_args.add_argument('--sqrt_mass', action='store_true',
                           help='Use a square-root transform on the probability mass axis '
                                '(default: no mass axis transform)')
    view_args.add_argument('--result', type=int, default=1,
                           help='Number of result to plot (used when there are multiple '
                                'possible results for the pair, default: DEFAULT)')
    view_args.add_argument('--vertical_colour', type=str, default='#4859a0',
                           help='Hex colour for vertical inheritance')
    view_args.add_argument('--horizontal_colour', type=str, default='#c47e7e',
                           help='Hex colour for horizontal inheritance')
    view_args.add_argument('--ambiguous_colour', type=str, default='#c9c9c9',
                           help='Hex colour for ambiguous inheritance')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='Verticall v' + __version__,
                            help="Show program's version number and exit")


def pairwise_and_view_settings(group):
    """
    The pairwise and view subcommands share a lot of settings in common, defined in this function.
    """
    setting_args = group.add_argument_group('Settings')
    setting_args.add_argument('--window_count', type=int, default=50000,
                              help='Aim to have at least this many comparison windows between '
                                   'assemblies')
    setting_args.add_argument('--ignore_indels', action='store_true',
                              help='Only use mismatches to determine distance (default: use '
                                   'both mismatches and gap-compressed indels)')
    setting_args.add_argument('--smoothing_factor', type=float, default=0.8,
                              help='Degree to which the distance distribution is smoothed')
    setting_args.add_argument('--secondary', type=float, default=0.7,
                              help='Peaks with a mass of at least this fraction of the most '
                                   'massive peak will be used to produce secondary distances')
    setting_args.add_argument('--verbose', action='store_true',
                              help='Output more detail to stderr for debugging (default: only '
                                   'output basic information)')

    alignment_args = group.add_argument_group('Alignment')
    # TODO: explore different indexing options (e.g. -k and -w) to see how they affect the results.
    alignment_args.add_argument('--index_options', type=str, default='-k15 -w10',
                                help='Minimap2 options for assembly indexing')
    # TODO: explore different alignment options (e.g. the things set by -x asm20).
    alignment_args.add_argument('--align_options', type=str, default='-x asm20',
                                help='Minimap2 options for assembly-to-assembly alignment')
    alignment_args.add_argument('--allowed_overlap', type=int, default=100,
                                help='Allow this much overlap between alignments')


def matrix_subparser(subparsers):
    group = subparsers.add_parser('matrix', description='produce a PHYLIP distance matrix',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('-i', '--in_file', type=pathlib.Path, required=True,
                               help='Filename of tsv created by vertical pairwise')
    required_args.add_argument('-o', '--out_file', type=pathlib.Path, required=True,
                               help='Filename for PHYLIP matrix output')

    settings_args = group.add_argument_group('Settings')
    settings_args.add_argument('--distance_type', type=str, default='median_vertical',
                               choices=['mean', 'median', 'peak', 'mean_vertical',
                                        'median_vertical'],
                               help='Which distance to use in matrix (choose from mean, median, '
                                    'peak, vert_mean and vert_median, default: DEFAULT')
    settings_args.add_argument('--asymmetrical', action='store_true',
                               help='Do not average pairs to make symmetrical matrices (default: '
                                    'make matrices symmetrical)')
    settings_args.add_argument('--no_jukes_cantor', action='store_true',
                               help='Do not apply Jukes-Cantor correction (default: apply '
                                    'Jukes-Cantor correction)')
    settings_args.add_argument('--multi', type=str, default='first',
                               choices=['first', 'low', 'high'],
                               help='Behaviour when there are multiple results for a sample pair')
    settings_args.add_argument('--names', type=str,
                               help='Samples names to include in matrix (comma-delimited, '
                                    'default: include all samples)')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='Verticall v' + __version__,
                            help="Show program's version number and exit")


def mask_subparser(subparsers):
    group = subparsers.add_parser('mask', description='mask horizontal regions from a SNV matrix',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('-i', '--in_file', type=pathlib.Path, required=True,
                               help='Filename of tsv created by vertical pairwise')
    required_args.add_argument('-m', '--in_matrix', type=pathlib.Path, required=True,
                               help='Filename for SNV matrix to be masked')
    required_args.add_argument('-o', '--out_file', type=pathlib.Path, required=True,
                               help='Filename for masked SNV matrix output')

    settings_args = group.add_argument_group('Settings')
    settings_args.add_argument('--multi', type=str, default='first',
                               choices=['first', 'low', 'high'],
                               help='Behaviour when there are multiple results for a sample pair')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='Verticall v' + __version__,
                            help="Show program's version number and exit")


def summary_subparser(subparsers):
    group = subparsers.add_parser('summary', description='summarise regions for one assembly',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('-i', '--in_file', type=pathlib.Path, required=True,
                               help='Filename of tsv created by vertical pairwise')
    required_args.add_argument('-a', '--assembly', type=pathlib.Path, required=True,
                               help='Filename for assembly to be summarised')

    settings_args = group.add_argument_group('Settings')
    settings_args.add_argument('--all', action='store_true',
                               help='Output one line for all assembly positions (default: omit '
                                    'redundant adjacent lines)')
    settings_args.add_argument('--plot', action='store_true',
                               help='Instead of outputting a table, display an interactive plot')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='Verticall v' + __version__,
                            help="Show program's version number and exit")


def repair_subparser(subparsers):
    group = subparsers.add_parser('repair', description='repair assembly for use in Verticall',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('-i', '--in_file', type=pathlib.Path, required=True,
                               help='Filename of assembly in need of repair')
    required_args.add_argument('-o', '--out_file', type=pathlib.Path, required=True,
                               help='Filename for repaired assembly output (if the same as -i, '
                                    'the input file will be overwritten)')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version', version='Verticall v' + __version__,
                            help="Show program's version number and exit")


def check_pairwise_args(args):
    pass


def check_view_args(args):
    samples = args.names.split(',')
    if len(samples) != 2:
        sys.exit('Error: two sample names (comma-delimited) must be supplied to --names')
    if args.result < 1:
        sys.exit('Error: --result must be a positive integer')
    if not check_hex_colour(args.vertical_colour):
        sys.exit('Error: --vertical_colour must be a valid hex colour code (with leading #)')
    if not check_hex_colour(args.horizontal_colour):
        sys.exit('Error: --horizontal_colour must be a valid hex colour code (with leading #)')
    if not check_hex_colour(args.ambiguous_colour):
        sys.exit('Error: --ambiguous_colour must be a valid hex colour code (with leading #)')


def check_matrix_args(args):
    pass


def check_mask_args(args):
    pass


def check_summary_args(args):
    pass


def check_repair_args(args):
    pass


if __name__ == '__main__':
    main()
