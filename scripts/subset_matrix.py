#!/usr/bin/env python3
"""
This is a script which takes two arguments:
* a filename of a PHYLIP matrix
* a comma-delimited list of sample names

And it outputs (to stdout) a PHYLIP matrix containing only the specified sample names.

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


def get_arguments():
    parser = argparse.ArgumentParser(description='Subset a PHYLIP matrix')

    parser.add_argument('phylip', type=str,
                        help='Filename of full PHYLIP matrix')
    parser.add_argument('samples', type=str,
                        help='Comma-delimited list of sample names')
    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    matrix = load_matrix(args.phylip)
    subset_sample_names = sorted(args.samples.split(','))
    print(len(subset_sample_names))
    for name_a in subset_sample_names:
        print(name_a, end='')
        for name_b in subset_sample_names:
            distance = matrix[(name_a, name_b)]
            print(f'\t{distance}', end='')
        print()


def load_matrix(phylip_filename):
    """
    Loads a PHYLIP distance matrix into a dictionary with key=pair value=distance.
    """
    sample_names, all_distances = [], []
    with open(phylip_filename, 'rt') as phylip_file:
        for i, line in enumerate(phylip_file):
            if i == 0:
                continue
            parts = line.strip().split('\t')
            sample_names.append(parts[0])
            all_distances.append(parts[1:])
    matrix = {}
    for name_a, distances in zip(sample_names, all_distances):
        assert len(distances) == len(sample_names)
        for name_b, distance in zip(sample_names, distances):
            matrix[(name_a, name_b)] = distance
    return matrix


if __name__ == '__main__':
    main()
