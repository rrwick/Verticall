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
from multiprocessing import Pool
import numpy as np
import statistics
import sys

from .log import log, section_header, explanation


def distance(args):
    distances, aligned_fractions, sample_names = \
        load_distances(args.alignment_results, args.method, args.threads)
    add_self_distances(distances, aligned_fractions, sample_names)
    check_matrix_size(distances, sample_names, args.alignment_results)
    correct_distances(distances, aligned_fractions, sample_names, args.correction)
    if not args.asymmetrical:
        make_symmetrical(distances, sample_names)
    output_phylip_matrix(distances, sample_names)


def load_distances(alignment_results, method, threads):
    section_header('Finding distances from distributions')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')
    distances, aligned_fractions, sample_names = {}, {}, set()

    arg_list = []
    with open(alignment_results, 'rt') as results:
        for line in results:
            if not line.startswith('#'):
                arg_list.append((line, method))
    count = len(arg_list)

    # If only using a single thread, do the alignment in a simple loop (easier for debugging).
    if threads == 1:
        for a in arg_list:
            a_1, a_2, d, aligned_fraction = load_one_distance(a)
            sample_names.update([a_1, a_2])
            distances[(a_1, a_2)] = d
            aligned_fractions[(a_1, a_2)] = aligned_fraction
            log(f'({len(distances)}/{count}) {a_1} vs {a_2}: {d:.9f}')

    # If using multiple threads, do the alignments in a process pool.
    else:
        with Pool(processes=threads) as pool:
            for a_1, a_2, d, aligned_fraction in pool.imap(load_one_distance, arg_list):
                sample_names.update([a_1, a_2])
                distances[(a_1, a_2)] = d
                aligned_fractions[(a_1, a_2)] = aligned_fraction
                log(f'({len(distances)}/{count}) {a_1} vs {a_2}: {d:.9f}')

    log()
    return distances, aligned_fractions, sorted(sample_names)


def load_one_distance(all_args):
    line, method = all_args
    parts = line.strip().split('\t')
    assembly_1, assembly_2 = parts[0], parts[1]
    piece_size = int(parts[2])
    assert piece_size > 0
    aligned_fraction = float(parts[3])
    assert 0.0 <= aligned_fraction <= 1.0
    masses = [float(p) for p in parts[4:]]
    d = get_distance(masses, piece_size, method)
    return assembly_1, assembly_2, d, aligned_fraction


def get_distance(masses, piece_size, method):
    if method == 'mean':
        d = get_mean(masses)
    elif method == 'median':
        d = get_interpolated_median(masses)
    elif method == 'mode':
        d = get_mode(masses)
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


def smooth_distribution(masses):
    """
    Smooths the distribution by redistributing mass between neighbouring points. Equal amounts of
    mass are moved up and down the distribution, so this smoothing doesn't change the mean.

    More smoothing is applied to higher masses, so the very low end of the distribution should
    remain relatively unchanged.
    """
    masses.append(0.0)
    changes = [0.0] * (len(masses))
    for i, m in enumerate(masses):
        if i == 0 or i == len(masses)-1:
            continue
        share_fraction = (2 ** (-100/i)) / 2
        share_amount = masses[i] * share_fraction
        changes[i-1] += share_amount / 2.0
        changes[i] -= share_amount
        changes[i+1] += share_amount / 2.0
    return [m+c for m, c in zip(masses, changes)]


def get_peak_distance(masses, max_tries=250):
    """
    Starts with the median and climb upwards, smoothing when a peak is reached. If a lot of
    smoothing doesn't change the peak, the process stops.

    The final returned value is interpolated from the peak and its neighbours.
    """
    median = get_median(masses)
    peak = climb_to_peak(masses, median)
    tries = 0
    while True:
        # If our peak is at zero, then no further smoothing is needed, so we can stop to save time.
        if peak == 0:
            break

        peak_before_smoothing = peak
        masses = smooth_distribution(masses)
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
