"""
This module contains code related to making a distance matrix and applying distance corrections.

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
import math

from .log import log, section_header, explanation


def save_all_matrices(args, sample_names, all_distances):
    section_header('Saving distance matrices')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')
    save_matrix(args, sample_names, all_distances, 'mean')
    save_matrix(args, sample_names, all_distances, 'median')
    save_matrix(args, sample_names, all_distances, 'mean_vertical')
    save_matrix(args, sample_names, all_distances, 'median_vertical')
    log()


def save_matrix(args, sample_names, all_distances, distance_type):
    distances, aligned_fractions = {}, {}
    for a in sample_names:
        for b in sample_names:
            if a == b:
                distances[(a, b)] = 0.0
                aligned_fractions[(a, b)] = 1.0
            else:
                distances[(a, b)] = all_distances[(a, b)][distance_type]
                aligned_fractions[(a, b)] = all_distances[(a, b)]['aligned_frac']

    correct_distances(distances, aligned_fractions, sample_names, args.correction)
    if not args.asymmetrical:
        make_symmetrical(distances, sample_names)

    output_filename = args.out_dir / (distance_type + '.phylip')
    log(f'Saving {output_filename}')
    with open(output_filename, 'wt') as f:
        f.write(str(len(sample_names)))
        f.write('\n')
        for a in sample_names:
            f.write(a)
            for b in sample_names:
                f.write(f'\t{distances[(a, b)]:.8f}')
            f.write('\n')


def correct_distances(distances, aligned_fractions, sample_names, correction):
    if 'jukescantor' in correction:
        for a in sample_names:
            for b in sample_names:
                distances[(a, b)] = jukes_cantor(distances[(a, b)])
    if 'alignedfrac' in correction:
        for a in sample_names:
            for b in sample_names:
                distances[(a, b)] = distances[(a, b)] / aligned_fractions[(a, b)]


def jukes_cantor(d):
    """
    https://www.desmos.com/calculator/okovk3dipx
    """
    if d == 0.0:
        return 0.0
    if d >= 0.75:
        return 25.0
    return -0.75 * math.log(1.0 - 1.3333333333333 * d)


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
