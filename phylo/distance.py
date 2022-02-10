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
from scipy.stats import nbinom
import statistics
import sys

from .gamma import fit_gamma_to_distribution
from .misc import get_mean
from .negbin import fit_negbin_to_distribution
from .log import log


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
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            assembly_1, assembly_2 = parts[0], parts[1]
            count = len(distances) + 1
            log(f'{count}: {assembly_1} vs {assembly_2}: ', end='')
            piece_size = int(parts[2])
            masses = [float(p) for p in parts[4:]]
            sample_names.update([assembly_1, assembly_2])
            d = get_distance(masses, piece_size, method)
            log(d)
            distances[(assembly_1, assembly_2)] = d
    log()
    return distances, sorted(sample_names)


def get_distance(masses, piece_size, method):
    if method == 'mean':
        d = get_mean(masses)
    elif method == 'median':
        d = get_median(masses)
    elif method == 'median_int':
        d = get_median_int(masses)
    elif method == 'mode':
        d = get_mode(masses)
    elif method == 'tightest_half':
        low, high = get_tightest_half(masses)
        d = statistics.mean([low, high])
    elif method == 'gamma':
        d = get_gamma_fit_distance(masses)
    elif method == 'negbin':
        d = get_negbin_fit_distance(masses)
    else:
        assert False
    return d / piece_size


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


def get_median_int(masses):
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


def get_tightest_half(masses):
    """
    Returns the low and high bounds which capture half (or more) of the total mass with the
    smallest difference between low and high.
    Results are zero-based and inclusive, e.g. (1, 2) means that the 2nd and 3rd positions together
    contain at least half of the mass.
    """
    half_total_mass = sum(masses) / 2.0
    best_low, best_high, best_range_size, best_total_in_range = None, None, None, None
    for low in range(0, get_median(masses) + 1):
        for high in range(low, len(masses)):
            total_in_range = sum(m for m in masses[low:high+1])
            if total_in_range >= half_total_mass:
                break
        range_size = high - low
        if ((best_range_size is None) or (range_size < best_range_size) or
                (range_size == best_range_size and total_in_range > best_total_in_range)):
            best_range_size = range_size
            best_total_in_range = total_in_range
            best_low = low
            best_high = high
    return best_low, best_high


def get_gamma_fit_distance(masses):
    """
    Fits a gamma distribution to the empirical distribution and returns the mean.
    """
    shape, scale, _ = fit_gamma_to_distribution(masses)
    gamma_mean = shape * scale
    return gamma_mean


def get_negbin_fit_distance(masses):
    """
    Fits a negative binomial distribution to the empirical distribution and returns the mean.
    """
    n, p, _ = fit_negbin_to_distribution(masses)
    return nbinom.mean(n, p)


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
