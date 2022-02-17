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
    elif method == 'top_quarter':
        d = get_top_quarter_median_distance(masses)
    elif method == 'peak':
        d, _ = get_peak_distance(masses)
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


def get_top_fraction(masses, fraction):
    """
    Returns low and high bounds which capture half (or more) of the total mass. The range starts
    with the median and climbs the distribution (shifting left or right to get a larger mass) and
    greedily expands the range. The range is returned in a Pythonic manner (0-based, exclusive-end).


    Since these results can be used for a mean/median, we don't want to return a result of (0, 1),
    i.e. a single value at the zero point of the distribution (will be common with very closely
    related genomes), as this will lead to a distance of zero. So in this situation we extend the
    high end of the distribution as long as it drops.
    """
    target_mass = sum(masses) * fraction
    median = get_median(masses)

    low, high = median, median+1
    while sum(masses[low:high]) < target_mass:

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
    low, high = get_top_fraction(masses, 0.5)
    masked_masses = [m if low <= i < high else 0.0 for i, m in enumerate(masses)]
    return get_mean(masked_masses)


def get_top_half_median_distance(masses):
    low, high = get_top_fraction(masses, 0.5)
    masked_masses = [m if low <= i < high else 0.0 for i, m in enumerate(masses)]
    return get_interpolated_median(masked_masses)


def get_top_quarter_mean_distance(masses):
    low, high = get_top_fraction(masses, 0.25)
    masked_masses = [m if low <= i < high else 0.0 for i, m in enumerate(masses)]
    return get_mean(masked_masses)


def get_top_quarter_median_distance(masses):
    low, high = get_top_fraction(masses, 0.25)
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


def smooth_distribution(masses, iterations):
    """
    Smooths the distribution by redistributing mass between neighbouring points. Equal amounts of
    mass are moved up and down the distribution, so this smoothing doesn't change the mean. More
    smoothing is applied to higher masses, so the very low end of the distribution should remain
    relatively unchanged.
    """
    for i in range(iterations):
        masses = smooth_distribution_one_iteration(masses, 0.5)
    return masses


def smooth_distribution_one_iteration(masses, max_share):
    masses.append(0.0)
    changes = [0.0] * (len(masses))
    mass_count = max(len(masses), 10)
    for i, m in enumerate(masses):
        if i == 0 or i == len(masses)-1:
            continue
        share_fraction = max_share * i / mass_count
        share_amount = masses[i] * share_fraction
        changes[i-1] += share_amount / 2.0
        changes[i] -= share_amount
        changes[i+1] += share_amount / 2.0
    return [m+c for m, c in zip(masses, changes)]


def get_peak_distance(masses, max_tries=1000):
    """
    Starts with the median and climb upwards, smoothing when a peak is reached. If a lot of
    smoothing doesn't change the peak, the process stops.

    The final returned value is interpolated from the peak and its neighbours.
    """
    median = get_median(masses)
    peak = climb_to_peak(masses, median)
    tries = 0
    while True:
        # If our peak is the global maximum, then there's no need to look further.
        if masses[peak] == max(masses):
            break

        peak_before_smoothing = peak
        masses = smooth_distribution(masses, 1)
        peak = climb_to_peak(masses, peak)
        if peak == peak_before_smoothing:  # if the smoothing didn't change the peak
            tries += 1
        else:  # if we got a new peak after smoothing
            tries = 0
        if tries == max_tries:
            break
    adjustment = interpolate(masses[peak-1] if peak > 0 else 0.0,
                             masses[peak],
                             masses[peak+1] if peak < len(masses)-1 else 0.0)
    return peak + adjustment, masses


def climb_to_peak(masses, starting_point):
    peak = starting_point
    while True:
        lower_mass = masses[peak-1] if peak > 0 else float('-inf')
        higher_mass = masses[peak+1] if peak < len(masses)-1 else float('-inf')
        if lower_mass >= masses[peak] and lower_mass > higher_mass:
            peak -= 1
        elif higher_mass > masses[peak] and higher_mass > lower_mass:
            peak += 1
        else:
            break
    return peak


def interpolate(low, peak, high):
    """
    This function takes three masses as input: the mass below the peak, the mass at the peak and
    the mass above the peak (i.e. it assumes that the below/above masses are no larger than the
    peak mass). It then returns an interpolation adjustment that varies from -0.5 to 0.5, similar
    to how interpolated medians work.
    """
    min_mass = min(low, peak, high)
    low -= min_mass
    peak -= min_mass
    high -= min_mass
    try:
        return (high - low) / (2.0 * peak)
    except ZeroDivisionError:
        return 0.0
