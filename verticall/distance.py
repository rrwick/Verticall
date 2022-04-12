"""
This module contains code related to making a distance distribution using the assembly-vs-assembly
alignments.

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
import math
import numpy as np
import statistics


def get_distribution(args, alignments):
    """
    Uses the alignments to build a distance distribution.
    """
    all_cigars = [a.simplified_cigar for a in alignments]
    window_size, window_step = choose_window_size_and_step(all_cigars, args.window_count)
    for a in alignments:
        a.set_up_sliding_windows(window_size, window_step)

    distances = []
    for a in alignments:
        distances += a.window_differences

    if len(distances) == 0:
        log_text = [f'V  no distances sampled']
        return None, window_size, len(distances), None, None, log_text

    distance_counts = collections.Counter(distances)
    masses = [0 if distance_counts[i] == 0 else distance_counts[i] / len(distances)
              for i in range(max(distances) + 1)]
    mean_distance = get_distance(masses, window_size, 'mean')
    median_distance = get_distance(masses, window_size, 'median')

    log_text = [f'V  distances sampled from sliding windows:',
                f'V    window size: {window_size} bp',
                f'V    window count: {len(distances)}',
                f'V    mean distance:   {mean_distance:.9f}',
                f'V    median distance: {median_distance:.9f}']

    return masses, window_size, len(distances), mean_distance, median_distance, log_text


def get_vertical_horizontal_distributions(alignments):
    """
    Returns two mass distributions, one for the vertical windows and another for the horizontal
    windows. Assumes that the alignments have already been painted.
    """
    vertical_distances, horizontal_distances = [], []
    for a in alignments:
        vertical_distances += a.get_all_vertical_distances()
        horizontal_distances += a.get_all_horizontal_distances()

    max_vertical_distance = max(vertical_distances) if vertical_distances else 0
    max_horizontal_distance = max(horizontal_distances) if horizontal_distances else 0
    max_distance = max(max_vertical_distance, max_horizontal_distance)

    total_length = len(vertical_distances) + len(horizontal_distances)
    vertical_counts = collections.Counter(vertical_distances)
    horizontal_counts = collections.Counter(horizontal_distances)
    vertical_masses = [0 if vertical_counts[i] == 0 else vertical_counts[i] / total_length
                       for i in range(max_distance + 1)]
    horizontal_masses = [0 if horizontal_counts[i] == 0 else horizontal_counts[i] / total_length
                         for i in range(max_distance + 1)]
    return vertical_masses, horizontal_masses


def choose_window_size_and_step(cigars, target_window_count):
    """
    This function chooses an appropriate window size and step for the given CIGARs. It tries to
    balance larger windows, which give higher-resolution identity samples, especially with
    closely-related assemblies, and smaller windows, which allow for more identity samples.
    """
    window_step = 1000
    while window_step > 1:
        window_size = window_step * 100
        if get_sliding_window_count(cigars, window_size, window_step) > target_window_count:
            return window_size, window_step
        window_step -= 1
    return window_step * 100, window_step


def get_sliding_window_count(cigars, window_size, window_step):
    """
    For a given window size, window step and set of CIGARs, this function returns how many windows
    there will be in total.
    """
    count = 0
    for cigar in cigars:
        cigar_len = len(cigar)
        if cigar_len < window_size:
            continue
        cigar_len -= window_size
        count += 1
        count += cigar_len // window_step
    return count


def get_distance(masses, piece_size, method):
    if method == 'mean':
        d = get_mean(masses)
    elif method == 'median':
        d = get_interpolated_median(masses)
    elif method == 'mode':
        d = get_mode(masses)
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


def get_peak_distance(masses, window_size, secondary_ratio):
    """
    Takes smoothed masses as input and finds the peak with the most mass. If there are multiple
    peaks closest to the most massive, then other peaks will also be outputted as secondary peaks.

    Returns:
    * a string listing all mass peaks (comma-delimited)
    * a list of results, each item a tuple of:
       * mass of the peak
       * result level (primary or secondary)
       * interpolated peak distance
       * thresholds for alignment painting
    * a list of log text
    """
    if masses is None:
        return None, [(None, None, None, None)], []

    log_text = ['V  mass peaks:']
    peaks_with_mass = [(get_peak_total_mass(masses, p), p) for p in find_peaks(masses)]
    largest_mass, most_massive_peak = sorted(peaks_with_mass)[-1]
    secondary_threshold = secondary_ratio * largest_mass

    mass_peaks, used_peaks = [], []
    for mass, peak in peaks_with_mass:
        mass_peak = f'{peak / window_size:.9f}'
        mass_peaks.append(mass_peak)
        if peak == most_massive_peak:
            result_level = 'primary'
            used_peaks.append((mass, peak, result_level))
        elif mass >= secondary_threshold:
            result_level = 'secondary'
            used_peaks.append((mass, peak, result_level))
        else:
            result_level = ''
        note = '' if not result_level else f' <- {result_level}'
        log_text.append(f'V    {mass_peak} ({100.0 * mass:.1f}%){note}')
    mass_peaks = ','.join(mass_peaks)

    results = []
    for mass, peak_count, result_level in used_peaks:
        thresholds = get_thresholds(masses, peak_count)
        mass_below = masses[peak_count-1] if peak_count > 0 else 0.0
        mass_at = masses[peak_count]
        mass_above = masses[peak_count+1] if peak_count < len(masses)-1 else 0.0
        peak_distance = (peak_count + interpolate(mass_below, mass_at, mass_above)) / window_size
        log_text.append(f'V    interpolated peak distance: {peak_distance:.9f} <- {result_level}')
        results.append((mass, result_level, peak_distance, thresholds))

    results = sorted(results, reverse=True)  # sort results from most to least massive
    return mass_peaks, results, log_text


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
    try:
        return (high - low) / (2.0 * (peak - min(low, peak, high)))
    except ZeroDivisionError:
        return 0.0


def find_peaks(masses):
    """
    Given a mass distribution, this returns a list of all peaks indices.
    """
    peaks = []
    for i, m in enumerate(masses):

        # If the point is not greater than the left-neighbour, move on.
        if not (True if i == 0 else m > masses[i-1]):
            continue

        # If the point is also greater than the right-neighbour, it's definitely a peak.
        if True if i == len(masses)-1 else m > masses[i+1]:
            peaks.append(i)
            continue

        # If the point is equal to its right-neighbour, then we need to check whether there is a
        # multi-point peak. If so, we add the middle position of the multi-point peak (rounded
        # down).
        if m == masses[i+1]:
            j = i
            while j < len(masses) and masses[j] == m:
                j += 1
            if j == len(masses) or m > masses[j]:
                peaks.append(int((i+j-1)/2))

    return peaks


def get_peak_total_mass(masses, peak):
    """
    Given a mass distribution and peak index, this returns the total mass of the peak by extending
    in both directions until the masses rise.
    """
    total = masses[peak]

    # Add masses on the left side of the peak.
    previous_mass = masses[peak]
    low = peak-1
    while low >= 0 and masses[low] <= previous_mass:
        total += masses[low]
        previous_mass = masses[low]
        low -= 1

    # Add masses on the right side of the peak.
    previous_mass = masses[peak]
    high = peak+1
    while high < len(masses) and masses[high] <= previous_mass:
        total += masses[high]
        previous_mass = masses[high]
        high += 1

    return total


def get_thresholds(masses, peak):
    low, very_low = get_low_thresholds(masses, peak)
    high, very_high = get_high_thresholds(masses, peak)
    return {'low': low, 'very_low': very_low, 'high': high, 'very_high': very_high}


def get_low_thresholds(masses, peak):
    minimum = find_local_minimum_to_left(masses, peak)
    if minimum is None:
        return None, None
    low_peak = find_local_maximum_to_left(masses, minimum)
    if low_peak is None:
        low_peak = 0
    low_1 = (peak + minimum) / 2
    low_2 = (minimum + low_peak) / 2
    return low_1, low_2


def get_high_thresholds(masses, peak):
    minimum = find_local_minimum_to_right(masses, peak)
    if minimum is None:
        return None, None
    high_peak = find_local_maximum_to_right(masses, minimum)
    if high_peak is None:
        high_peak = len(masses)-1
    high_1 = (peak + minimum) / 2
    high_2 = (minimum + high_peak) / 2
    return high_1, high_2


def find_local_minimum_to_right(masses, i):
    """
    Starting at a given index, this function looks for a local minimum to the right. If one exists,
    its index is returned. If one does not exist (e.g. the masses decrease all the way to the end),
    then None is returned.
    """
    if i == len(masses) - 1:
        return None
    while masses[i+1] <= masses[i]:
        i += 1
        if i == len(masses)-1:
            return None
    return i


def find_local_minimum_to_left(masses, i):
    """
    Starting at a given index, this function looks for a local minimum to the right. If one exists,
    its index is returned. If one does not exist (e.g. the masses decrease all the way to the end),
    then None is returned.
    """
    if i == 0:
        return None
    while masses[i-1] <= masses[i]:
        i -= 1
        if i == 0:
            return None
    return i


def find_local_maximum_to_right(masses, i):
    """
    Starting at a given index, this function looks for a local maximum to the right. If one exists,
    its index is returned. If one does not exist (e.g. the masses increase all the way to the end),
    then None is returned.
    """
    if i == len(masses) - 1:
        return None
    while masses[i+1] > masses[i]:
        i += 1
        if i == len(masses)-1:
            return None
    return i


def find_local_maximum_to_left(masses, i):
    """
    Starting at a given index, this function looks for a local maximum to the right. If one exists,
    its index is returned. If one does not exist (e.g. the masses increase all the way to the end),
    then None is returned.
    """
    if i == 0:
        return None
    while masses[i-1] > masses[i]:
        i -= 1
        if i == 0:
            return None
    return i


def smooth_distribution(masses, smoothing_factor):
    if masses is None:
        return None

    smoothed = []
    for i, _ in enumerate(masses):
        kernel_width = i ** smoothing_factor
        smoothed.append(get_smoothed_mass(masses, i, kernel_width))

    # Normalise to sum to one.
    total = sum(smoothed)
    return [s/total for s in smoothed]


def get_smoothed_mass(masses, i, kernel_width):
    low_i = max(math.floor(i - kernel_width), 0)
    high_i = math.ceil(i + kernel_width)
    masses_to_average, weights = [], []
    for j in range(low_i, high_i+1):
        try:
            masses_to_average.append(masses[j])
        except IndexError:
            masses_to_average.append(0.0)
        weights.append(get_epanechnikov_weight(kernel_width, j-i))
    return np.average(masses_to_average, weights=weights)


def get_epanechnikov_weight(kernel_width, offset):
    """
    Returns the weight of an Epanechnikov kernel. Since these weights are just relative, we don't
    bother normalising the kernel to an area of 1.
    """
    if kernel_width == 0.0:
        return 1.0 if offset == 0.0 else 0.0
    weight = 1.0 - ((offset / kernel_width) ** 2)
    return max(0.0, weight)
