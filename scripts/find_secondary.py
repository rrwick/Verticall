#!/usr/bin/env python3
"""
This is a script which takes two or three arguments:
* a filename of a Verticall TSV file
* a comma-delimited list of sample names
* (optional) another comma-delimited list of sample names

If given one list of sample names, it outputs (to stdout) the lines of the TSV where both of the
samples were in the given list and a secondary result exists for that pair.

If given two lists of sample names, it outputs (to stdout) the lines of the TSV where both of the
samples were in the given list and a secondary result exists for that pair.

This script can be useful for making manual edits to the TSV file in the case of secondary results.

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
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Subset a PHYLIP matrix')

    parser.add_argument('tsv', type=str,
                        help='Filename of Verticall tsv file')
    parser.add_argument('samples', type=str, nargs='+',
                        help='Comma-delimited list of sample names')
    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    if len(args.samples) == 1:
        one_list(args.tsv, args.samples[0])
    elif len(args.samples) == 2:
        two_lists(args.tsv, args.samples[0], args.samples[1])
    else:
        sys.exit('Error: you must supply either one or two comma-delimited lists of sample names')


def one_list(tsv_filename, samples_str):
    samples = set(samples_str.split(','))
    level_index = None
    pairs_to_print = set()
    with open(tsv_filename, 'rt') as tsv:
        for i, line in enumerate(tsv):
            parts = line.strip().split('\t')
            if i == 0:
                print(line.strip())
                level_index = parts.index('result_level')
            elif parts[0] in samples and parts[1] in samples and parts[level_index] == 'secondary':
                pairs_to_print.add((parts[0], parts[1]))
    with open(tsv_filename, 'rt') as tsv:
        for i, line in enumerate(tsv):
            parts = line.strip().split('\t')
            if i == 0:
                continue
            if (parts[0], parts[1]) in pairs_to_print:
                print(line.strip())


def two_lists(tsv_filename, samples_str_1, samples_str_2):
    samples_1 = set(samples_str_1.split(','))
    samples_2 = set(samples_str_2.split(','))
    level_index = None
    pairs_to_print = set()
    with open(tsv_filename, 'rt') as tsv:
        for i, line in enumerate(tsv):
            parts = line.strip().split('\t')
            if i == 0:
                print(line.strip())
                level_index = parts.index('result_level')
            elif ((parts[0] in samples_1 and parts[1] in samples_2)
                  or (parts[0] in samples_2 and parts[1] in samples_1))\
                    and parts[level_index] == 'secondary':
                pairs_to_print.add((parts[0], parts[1]))
    with open(tsv_filename, 'rt') as tsv:
        for i, line in enumerate(tsv):
            parts = line.strip().split('\t')
            if i == 0:
                continue
            if (parts[0], parts[1]) in pairs_to_print:
                print(line.strip())


if __name__ == '__main__':
    main()
