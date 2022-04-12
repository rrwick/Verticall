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
import sys

from .log import log, section_header, explanation, warning
from .misc import check_file_exists


def matrix(args):
    welcome_message()
    distances, sample_names = load_tsv_file(args.in_file, args.distance_type)
    if args.names is not None:
        sample_names = filter_names(sample_names, args.names)
    if not args.no_jukes_cantor:
        jukes_cantor_correction(distances, sample_names)
    if not args.asymmetrical:
        make_symmetrical(distances, sample_names)
    save_matrix(args.out_file, distances, sample_names)
    finished_message()


def welcome_message():
    section_header('Starting Verticall matrix')
    explanation('Verticall matrix extracts distances from the tab-delimited file made by Vertical '
                'pairwise, producing a PHYLIP distance matrix suitable for use in distance-based '
                'phylogeny algorithms (e.g. BioNJ or FastME).')
    # TODO: display the options used here


def finished_message():
    section_header('Finished!')
    explanation('You can now use a distance-based algorithm (e.g. BioNJ or FastME) to produce a '
                'phylogeny from the distance matrix.')


def load_tsv_file(filename, distance_type):
    # TODO: add --multi logic
    check_file_exists(filename)
    section_header('Loading distances')
    distances, sample_names = {}, set()
    column_index = None
    with open(filename, 'rt') as f:
        for i, line in enumerate(f):
            parts = line.strip('\n').split('\t')
            if i == 0:  # header line
                column_index = get_column_index(parts, distance_type, filename)
            else:
                assembly_a, assembly_b = parts[0], parts[1]
                try:
                    if parts[column_index] == '':
                        distance = None
                    else:
                        distance = float(parts[column_index])
                except ValueError:
                    sys.exit(f'Error: could not convert {parts[column_index]} to a number')
                sample_names.add(assembly_a)
                sample_names.add(assembly_b)
                distances[(assembly_a, assembly_b)] = distance
    sample_names = sorted(sample_names)
    log(f'{len(distances)} distances loaded for {len(sample_names)} assemblies')
    for sample_name in sample_names:
        distances[(sample_name, sample_name)] = 0.0
    check_for_missing_distances(distances, sample_names)
    log()
    return distances, sorted(sample_names)


def filter_names(all_names, specified_names):
    all_names = set(all_names)
    filtered_names = set()
    for name in specified_names.split(','):
        if name in all_names:
            filtered_names.add(name)
        else:
            sys.exit(f'Error: could not find sample {name}')
    return sorted(filtered_names)


def check_for_missing_distances(distances, sample_names):
    any_missing = False
    for sample_a in sample_names:
        for sample_b in sample_names:
            if (sample_a, sample_b) not in distances:
                any_missing = True
                distances[(sample_a, sample_b)] = None
    if any_missing:
        warning('some pairwise distances were not found - matrix will contain empty cells.')


def get_column_index(header_parts, distance_type, filename):
    if header_parts[0] != 'assembly_a':
        sys.exit(f'Error: first column in {filename} is not labelled "assembly_a" - is the file '
                 f'formatted correctly?')
    if header_parts[1] != 'assembly_b':
        sys.exit(f'Error: second column in {filename} is not labelled "assembly_b" - is the file '
                 f'formatted correctly?')
    target_header = distance_type + '_distance'
    for i, header in enumerate(header_parts):
        if target_header == header:
            return i
    sys.exit(f'Error: could not find {target_header} column in {filename}')


def save_matrix(filename, distances, sample_names):
    section_header('Saving matrix to file')
    log(f'{filename.resolve()}')

    distance_count, missing_distances = 0, False
    with open(filename, 'wt') as f:
        f.write(str(len(sample_names)))
        f.write('\n')
        for a in sample_names:
            f.write(a)
            for b in sample_names:
                distance = distances[(a, b)]
                if distance is None:
                    distance = ''
                    missing_distances = True
                else:
                    distance = f'{distance:.9f}'
                    distance_count += 1
                f.write(f'\t{distance}')
            f.write('\n')

    log(f'{len(sample_names)} samples, {distance_count} distances')
    log()
    if missing_distances:
        log(bold_red('WARNING: '
                     'one or more distances are missing resulting in an incomplete matrix'))
        log()


def jukes_cantor_correction(distances, sample_names):
    for a in sample_names:
        for b in sample_names:
            distances[(a, b)] = jukes_cantor(distances[(a, b)])


def jukes_cantor(d):
    """
    https://www.desmos.com/calculator/okovk3dipx
    """
    if d is None:
        return None
    if d == 0.0:
        return 0.0
    if d >= 0.75:
        return 25.0
    return -0.75 * math.log(1.0 - 1.3333333333333 * d)


def make_symmetrical(distances, sample_names):
    for a, b in itertools.combinations(sample_names, 2):
        d1 = distances[(a, b)]
        d2 = distances[(b, a)]
        if d1 is not None and d2 is not None:
            mean_distance = (d1 + d2) / 2.0
        elif d1 is None and d2 is None:
            mean_distance = None
        elif d1 is None:
            mean_distance = d2
        else:
            mean_distance = d1
        distances[(a, b)] = mean_distance
        distances[(b, a)] = mean_distance
