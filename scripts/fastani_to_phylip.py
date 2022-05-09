#!/usr/bin/env python3
"""
This is a script which takes as input one or more fastANI files and outputs (to stdout) a PHYLIP
distance matrix. It averages the two orders of each pair (e.g. A-vs-B and B-vs-A) so the matrix is
symmetrical.

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

import itertools
import pathlib
import sys


def main():
    sample_names = set()
    distances = {}
    for fastani_file in sys.argv[1:]
        with fastani_file as f:
            for line in f:
                parts = line.strip().split('\t')
                a_1 = pathlib.Path(parts[0]).name.split('.fasta')[0]
                a_2 = pathlib.Path(parts[1]).name.split('.fasta')[0]
                ani = float(parts[2])
                distance = 1.0 - (ani / 100.0)
                sample_names.add(a_1)
                sample_names.add(a_2)
                distances[(a_1, a_2)] = distance
    sample_names = sorted(sample_names)
    make_symmetrical(distances, sample_names)
    output_phylip_matrix(distances, sample_names)


def make_symmetrical(distances, sample_names):
    for a, b in itertools.combinations(sample_names, 2):
        d1 = distances[(a, b)]
        d2 = distances[(b, a)]
        mean_distance = (d1 + d2) / 2.0
        distances[(a, b)] = mean_distance
        distances[(b, a)] = mean_distance


def output_phylip_matrix(distances, sample_names):
    print(len(sample_names))
    for a in sample_names:
        print(a, end='')
        for b in sample_names:
            print(f'\t{distances[(a, b)]:.8f}', end='')
        print()


if __name__ == '__main__':
    main()
