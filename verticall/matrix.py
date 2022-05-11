"""
This module contains code for the 'verticall matrix' subcommand.

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

import collections
import itertools
import math
import sys

from .log import log, section_header, explanation, warning
from .misc import check_file_exists
from .tsv import get_column_index


def matrix(args):
    welcome_message()
    distances, sample_names = load_tsv_file(args.in_file, args.distance_type)
    distances, sample_names = resolve_multi_distances(distances, sample_names, args.multi)
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
    explanation('Verticall matrix extracts distances from the tab-delimited file made by Verticall '
                'pairwise, producing a PHYLIP distance matrix suitable for use in distance-based '
                'phylogeny algorithms (e.g. BioNJ or FastME).')


def finished_message():
    section_header('Finished!')
    explanation('You can now use a distance-based algorithm (e.g. BIONJ or FastME) to produce a '
                'phylogeny from the distance matrix.')


def load_tsv_file(filename, distance_type):
    check_file_exists(filename)
    distances, sample_names = collections.defaultdict(list), set()
    column_index = None
    # excluded_samples = get_multi_result_samples(filename) if multi == 'exclude' else set()
    with open(filename, 'rt') as f:
        for i, line in enumerate(f):
            parts = line.strip('\n').split('\t')
            if i == 0:  # header line
                column_index = get_column_index(parts, distance_type + '_distance', filename)
            else:
                assembly_a, assembly_b = parts[0], parts[1]
                # if assembly_a in excluded_samples or assembly_b in excluded_samples:
                #     continue
                distance = get_distance_from_line_parts(parts, column_index)
                sample_names.add(assembly_a)
                sample_names.add(assembly_b)
                distances[(assembly_a, assembly_b)].append(distance)
    sample_names = sorted(sample_names)
    log(f'{len(distances)} distances loaded for {len(sample_names)} assemblies')
    for sample_name in sample_names:
        distances[(sample_name, sample_name)].append(0.0)
    log()
    return distances, sorted(sample_names)


def resolve_multi_distances(distances, sample_names, multi):
    section_header('Resolving multi-distance pairs')
    explanation('Some pairs may have more than one result in the TSV file, so this step reduces '
                'each pair to a single distance using the logic chosen by the --multi option.')

    multi_distance_pairs = set(p for p, d in distances.items() if len(d) > 1)
    log(f'Multi-distance pairs: {len(multi_distance_pairs)} / {len(distances)}')
    log()
    if len(multi_distance_pairs) == 0:
        log('No resolution required')
        return {p: d[0] for p, d in distances.items()}, sample_names
    elif multi == 'concordant':
        log('Resolving by greedily choosing the most tree-concordant distances')
        distances = choose_concordant_distances(distances, sample_names, multi_distance_pairs)
    elif multi == 'exclude':
        log('Resolving by excluding any samples in a multi-distance pair')
        distances, sample_names = exclude_multi_distances(distances, multi_distance_pairs)
    elif multi == 'first':
        log('Resolving using TSV file order (keeping the first distance for each pair)')
        distances = {p: d[0] for p, d in distances.items()}
    elif multi == 'low':
        log('Resolving to minimum (keeping the lowest distance for each pair)')
        distances =  {p: min(d) for p, d in distances.items()}
    elif multi == 'high':
        log('Resolving to maximum (keeping the highest distance for each pair)')
        distances = {p: max(d) for p, d in distances.items()}
    else:
        assert False
    log()
    return distances, sample_names


def choose_concordant_distances(distances, sample_names, multi_distance_pairs):
    round_num = 1
    while True:
        # Start with the first distance in each pair (like with --multi first) - this is the
        # baseline we will try to improve upon.
        first_distances = {p: d[0] for p, d in distances.items()}
        make_symmetrical(first_distances, sample_names)
        first_concordance = get_tree_concordance(first_distances)
        log(f'\nRound {round_num} starting concordance: {first_concordance}')

        # Now we go through each multi-distance pair and check the concordance using each possible
        # alternative distance, remembering the best result.
        best_concordance, best_pair_and_distance = None, None
        for pair in sorted(multi_distance_pairs):
            assert len(distances[pair]) > 1
            log(f'  {pair[0]} vs {pair[1]}:')
            for alt_distance in distances[pair][1:]:
                trial_distances = {p: d[0] for p, d in distances.items()}
                trial_distances[pair] = alt_distance
                concordance = get_tree_concordance(trial_distances)
                log(f'    {distances[pair][0]} -> {alt_distance}, concordance: {concordance}')
                if best_concordance is None or concordance < best_concordance:
                    best_concordance = concordance
                    best_pair_and_distance = (pair, alt_distance)
        log(f'  Round {round_num} best concordance:     {best_concordance}')

        # If our best result improves concordance, we cement that distance (remove alternatives).
        if best_concordance < first_concordance:
            best_pair, best_distance = best_pair_and_distance
            distances[best_pair] = [best_distance]
            multi_distance_pairs.remove(best_pair)
        else:
            log(f'\nNo improvement - multi-distance resolution finished')
            break
        if len(multi_distance_pairs) == 0:
            log(f'\nNo more multi-distance pairs remain')
            break

    return {p: d[0] for p, d in distances.items()}


def get_tree_concordance(distances):
    """
    This function takes in a distance matrix and returns a value which indicates how well the
    distances fit into a tree.
    """
    return 0.0  # TEMP


def exclude_multi_distances(distances, multi_distance_pairs):
    excluded_samples = set()
    for assembly_a, assembly_b in multi_distance_pairs:
        excluded_samples.add(assembly_a)
        excluded_samples.add(assembly_b)
    excluded_str = ', '.join(sorted(excluded_samples))
    log()
    log(f'The following samples will be excluded due to multiple results: {excluded_str}')
    new_distances, new_sample_names = {}, set()
    for pair, distance_list in distances.items():
        assembly_a, assembly_b = pair
        if assembly_a not in excluded_samples and assembly_b not in excluded_samples:
            assert len(distance_list) == 1
            new_distances[(assembly_a, assembly_b)] = distance_list[0]
            new_sample_names.add(assembly_a)
            new_sample_names.add(assembly_b)
    log()
    log(f'{len(new_sample_names)} samples ({len(new_distances)} distances) remain')
    return new_distances, sorted(new_sample_names)


def get_distance_from_line_parts(parts, column_index):
    try:
        distance = parts[column_index]
    except IndexError:
        sys.exit(f'Error: column {column_index+1} missing from tsv file')
    try:
        if distance == '':
            return None
        else:
            return float(distance)
    except ValueError:
        sys.exit(f'Error: could not convert {distance} to a number')


def multi_distance(existing_distance, new_distance, multi):
    """
    This code handles the case when a distance is seen a subsequent time for a single assembly pair.
    Given the existing distance (earlier in the tsv), a new distance (later in the tsv) and multi
    logic, it returns the distance that should be used in the matrix.
    """
    if multi == 'first':
        return existing_distance
    elif multi == 'low':
        return min(existing_distance, new_distance)
    elif multi == 'high':
        return max(existing_distance, new_distance)
    else:
        assert False


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
    log()
    if missing_distances:
        warning('one or more distances are missing resulting in an incomplete matrix')


def jukes_cantor_correction(distances, sample_names):
    """
    Applies Jukes-Cantor correction in-place to the entire distance matrix.
    """
    for a in sample_names:
        for b in sample_names:
            distances[(a, b)] = jukes_cantor(distances[(a, b)])


def jukes_cantor(d):
    """
    Applies Jukes-Cantor correction to a single distance, returning the corrected value.
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
    """
    Makes the distance matrix symmetrical, changing it in-place.
    """
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
