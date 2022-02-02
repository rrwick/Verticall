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

import itertools
import math
import numpy as np
import sys


def distance(args):
    distances, sample_names = load_distances(args.alignment_results, args.method)
    add_self_distances(distances, sample_names)
    check_matrix_size(distances, sample_names, args.alignment_results)
    correct_distances(distances, sample_names, args.correction)
    if not args.asymmetrical:
        make_symmetrical(distances, sample_names)
    output_phylip_matrix(distances, sample_names)


def load_distances(alignment_results, method):
    distances, sample_names = {}, set()
    with open(alignment_results, 'rt') as results:
        for line in results:
            parts = line.strip().split('\t')
            assembly_1, assembly_2 = parts[0], parts[1]
            piece_size = int(parts[2])
            masses = [float(p) for p in parts[4:]]
            sample_names.update([assembly_1, assembly_2])
            distances[(assembly_1, assembly_2)] = get_distance(masses, piece_size, method)
    return distances, sorted(sample_names)


def get_distance(masses, piece_size, method):
    if method == 'mean':
        return get_mean_distance(masses, piece_size)
    elif method == 'median':
        return get_median_distance(masses, piece_size)
    elif method == 'median_cont':
        return get_median_cont_distance(masses, piece_size)
    assert False


def get_mean_distance(masses, piece_size):
    distances = [i / piece_size for i in range(len(masses))]
    return np.average(distances, weights=masses)


def get_median_distance(masses, piece_size):
    return 0.0


def get_median_cont_distance(masses, piece_size):
    return 0.0


def add_self_distances(distances, sample_names):
    for n in sample_names:
        distances[(n, n)] = 0.0


def check_matrix_size(distances, sample_names, alignment_results):
    if len(distances) != len(sample_names)**2:
        sys.exit(f'\nError: incorrect number of records in {alignment_results}'
                 f' - rerun XXXXXXXXX align')


def correct_distances(distances, sample_names, correction):
    if correction == 'none':
        return
    elif correction == 'jukescantor':
        for a in sample_names:
            for b in sample_names:
                distances[(a, b)] = jukes_cantor(distances[(a, b)])
        return
    assert False


def jukes_cantor(d):
    if d == 0.0:
        return 0.0
    return -0.75 * math.log(1 - 1.33333333333333 * d)


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
