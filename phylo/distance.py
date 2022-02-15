"""
Copyright 2022 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/XXXXXXXXX

This file is part of XXXXXXXXX. XXXXXXXXX is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. XXXXXXXXX is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with XXXXXXXXX.
If not, see <https://www.gnu.org/licenses/>.
"""

import itertools
import math
import numpy as np
import statistics
import sys

from .log import log


def distance(args):
    distances, aligned_fractions, sample_names = load_distances(args.alignment_results, args.method)
    add_self_distances(distances, aligned_fractions, sample_names)
    check_matrix_size(distances, sample_names, args.alignment_results)
    correct_distances(distances, aligned_fractions, sample_names, args.correction)
    if not args.asymmetrical:
        make_symmetrical(distances, sample_names)
    output_phylip_matrix(distances, sample_names)


def load_distances(alignment_results, method):
    distances, aligned_fractions, sample_names = {}, {}, set()
    with open(alignment_results, 'rt') as results:
        for line in results:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            assembly_1, assembly_2 = parts[0], parts[1]
            count = len(distances) + 1
            piece_size = int(parts[2])
            assert piece_size > 0
            aligned_fraction = float(parts[3])
            assert 0.0 <= aligned_fraction <= 1.0
            masses = [float(p) for p in parts[4:]]
            sample_names.update([assembly_1, assembly_2])
            d = get_distance(masses, piece_size, method)
            log(f'{count}: {assembly_1} vs {assembly_2}: {d:.9f}')
            distances[(assembly_1, assembly_2)] = d
            aligned_fractions[(assembly_1, assembly_2)] = aligned_fraction
    log()
    return distances, aligned_fractions, sorted(sample_names)


def get_distance(masses, piece_size, method):
    if method == 'mean':
        d = get_mean(masses)
    elif method == 'median':
        d = get_interpolated_median(masses)
    elif method == 'mode':
        d = get_mode(masses)
    elif method == 'top_half':
        d = get_top_half_median_distance(masses)
    else:
        assert False
    return d / piece_size


def get_mean(masses):
    return np.average(range(len(masses)), weights=masses)


def get_median(masses):
    """
    Returns the median of the distance distribution. This median is not interpolated, i.e. it will
    be equal to one of the distances in the distribution.
    """
    half_total_mass = sum(masses) / 2.0
    total = 0.0
    for i, m in enumerate(masses):
        total += m
        if total >= half_total_mass:
            return i
    return 0


def get_interpolated_median(masses):
    """
    Returns the interpolated median of the distance distribution:
    https://en.wikipedia.org/wiki/Median#Interpolated_median
    https://aec.umich.edu/median.php
    """
    median = get_median(masses)
    below, equal, above = 0.0, 0.0, 0.0
    for i, m in enumerate(masses):
        if i < median:
            below += m
        elif i > median:
            above += m
        else:  # i == median
            equal += m
    if equal == 0.0:
        interpolated_median = median
    else:
        interpolated_median = median + ((above - below) / (2.0 * equal))
    return interpolated_median


def get_mode(masses):
    """
    Returns the mode of the distance distribution (the distance with the highest mass). If multiple
    distances tie for the highest mass (the distribution is multimodal), it returns the mean of
    those distances.
    """
    max_mass = max(masses)
    distances_with_max_mass = []
    for i, m in enumerate(masses):
        if m == max_mass:
            distances_with_max_mass.append(i)
    if len(distances_with_max_mass) == 1:
        return distances_with_max_mass[0]
    else:
        return statistics.mean(distances_with_max_mass)


def get_top_half(masses):
    """
    Returns low and high bounds which capture half (or more) of the total mass. The range starts
    with the median and climbs the distribution (shifting left or right to get a larger mass) and
    greedily expands the range. The range is returned in a Pythonic manner (0-based, exclusive-end).


    Since these results can be used for a mean/median, we don't want to return a result of (0, 1),
    i.e. a single value at the zero point of the distribution (will be common with very closely
    related genomes), as this will lead to a distance of zero. So in this situation we extend the
    high end of the distribution as long as it drops.
    """
    half_total_mass = sum(masses) / 2.0
    median = get_median(masses)

    low, high = median, median+1
    while sum(masses[low:high]) < half_total_mass:

        # If we've reached the limits on both ends (probably due to the min_samples values and a
        # small distribution), we're done.
        if low == 0 and high == len(masses):
            break

        # We first check to see if shifting the range up or down by one can increase the total.
        current_total = sum(masses[low:high])
        shift_down_total = sum(masses[low-1:high-1]) if low > 0 else float('-inf')
        shift_up_total = sum(masses[low+1:high+1]) if high < len(masses) else float('-inf')
        if shift_down_total > current_total and shift_down_total > shift_up_total:
            low -= 1
            high -= 1
            continue
        if shift_up_total > current_total and shift_up_total > shift_down_total:
            low += 1
            high += 1
            continue

        # If we got here, then shifting failed to increase the total, so we need to expand the
        # range. If we've reached the limit on either end, then we can only expand in one way.
        if low == 0:
            high += 1
            continue
        if high == len(masses):
            low -= 1
            continue

        # If we can potentially expand in either way, we need to decide which way to expand.
        if masses[high] > masses[low-1]:
            high += 1
            continue
        if masses[low-1] > masses[high]:
            low -= 1
            continue

        # If we got here, then the new low and new high are tied, so we expand in both directions.
        high += 1
        low -= 1

    # If the range starts at the bottom of the distribution, extend the top end until it either
    # hits zero or starts to rise.
    if low == 0:
        while high < len(masses) and 0.0 < masses[high] < masses[high-1]:
            high += 1

    return low, high


def get_top_half_mean_distance(masses):
    low, high = get_top_half(masses)
    masked_masses = [m if low <= i < high else 0.0 for i, m in enumerate(masses)]
    return get_mean(masked_masses)


def get_top_half_median_distance(masses):
    low, high = get_top_half(masses)
    masked_masses = [m if low <= i < high else 0.0 for i, m in enumerate(masses)]
    return get_interpolated_median(masked_masses)


def add_self_distances(distances, aligned_fractions, sample_names):
    for n in sample_names:
        distances[(n, n)] = 0.0
        aligned_fractions[(n, n)] = 1.0


def check_matrix_size(distances, sample_names, alignment_results):
    if len(distances) != len(sample_names)**2:
        sys.exit(f'\nError: incorrect number of records in {alignment_results}'
                 f' - rerun XXXXXXXXX align')


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
