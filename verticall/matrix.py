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

    aligned_fractions = get_matrix(sample_names, all_distances, 'aligned_fraction')
    mean_distances = get_matrix(sample_names, all_distances, 'mean')
    median_distances = get_matrix(sample_names, all_distances, 'median')
    peak_distances = get_matrix(sample_names, all_distances, 'peak')
    mean_vertical_distances = get_matrix(sample_names, all_distances, 'mean_vertical')
    median_vertical_distances = get_matrix(sample_names, all_distances, 'median_vertical')

    save_matrix(args, sample_names, aligned_fractions, 'aligned_fraction')
    save_matrix(args, sample_names, mean_distances, 'mean')
    save_matrix(args, sample_names, median_distances, 'median')
    save_matrix(args, sample_names, peak_distances, 'peak')
    save_matrix(args, sample_names, mean_vertical_distances, 'mean_vertical')
    save_matrix(args, sample_names, median_vertical_distances, 'median_vertical')

    correct_distances(median_vertical_distances, aligned_fractions, sample_names, args.correction)
    if not args.asymmetrical:
        make_symmetrical(median_vertical_distances, sample_names)

    save_matrix(args, sample_names, median_vertical_distances, 'final')
    log()


def get_matrix(sample_names, all_distances, distance_type):
    matrix = {}
    for a in sample_names:
        for b in sample_names:
            if a == b:
                if distance_type == 'aligned_fraction':
                    matrix[(a, b)] = 1.0
                else:
                    matrix[(a, b)] = 0.0
            else:
                matrix[(a, b)] = all_distances[(a, b)][distance_type]
    return matrix


def save_matrix(args, sample_names, matrix, file_prefix):
    output_filename = args.out_dir / (file_prefix + '.phylip')
    log(f'Saving {output_filename}')
    with open(output_filename, 'wt') as f:
        f.write(str(len(sample_names)))
        f.write('\n')
        for a in sample_names:
            f.write(a)
            for b in sample_names:
                f.write(f'\t{matrix[(a, b)]:.8f}')
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
