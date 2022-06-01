"""
This module contains some tests for Verticall. To run them, execute `pytest` from the root
Verticall directory.

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

import pytest

import verticall.distance


def test_get_mean():
    assert verticall.distance.get_mean([1.00, 0.00, 0.00, 0.00]) == pytest.approx(0.0)
    assert verticall.distance.get_mean([0.00, 1.00, 0.00, 0.00]) == pytest.approx(1.0)
    assert verticall.distance.get_mean([0.00, 0.00, 1.00, 0.00]) == pytest.approx(2.0)
    assert verticall.distance.get_mean([0.00, 0.00, 0.00, 1.00]) == pytest.approx(3.0)
    assert verticall.distance.get_mean([0.50, 0.50, 0.00, 0.00]) == pytest.approx(0.5)
    assert verticall.distance.get_mean([0.00, 0.50, 0.50, 0.00]) == pytest.approx(1.5)
    assert verticall.distance.get_mean([0.00, 0.00, 0.50, 0.50]) == pytest.approx(2.5)
    assert verticall.distance.get_mean([0.25, 0.25, 0.25, 0.25]) == pytest.approx(1.5)
    assert verticall.distance.get_mean([0.10, 0.20, 0.30, 0.40]) == pytest.approx(2.0)
    assert verticall.distance.get_mean([0.40, 0.30, 0.20, 0.10]) == pytest.approx(1.0)


def test_get_median():
    assert verticall.distance.get_median([1.0, 0.0, 0.0, 0.0]) == 0
    assert verticall.distance.get_median([0.0, 1.0, 0.0, 0.0]) == 1
    assert verticall.distance.get_median([0.0, 0.0, 1.0, 0.0]) == 2
    assert verticall.distance.get_median([0.0, 0.0, 0.0, 1.0]) == 3
    assert verticall.distance.get_median([0.6, 0.0, 0.0, 0.4]) == 0
    assert verticall.distance.get_median([0.4, 0.0, 0.0, 0.6]) == 3
    assert verticall.distance.get_median([0.1, 0.2, 0.3, 0.4]) == 2
    assert verticall.distance.get_median([0.4, 0.3, 0.2, 0.1]) == 1
    assert verticall.distance.get_median([]) == 0


def test_get_interpolated_median():
    assert verticall.distance.get_interpolated_median([1.0, 0.0, 0.0, 0.0]) == pytest.approx(0.0)
    assert verticall.distance.get_interpolated_median([0.0, 1.0, 0.0, 0.0]) == pytest.approx(1.0)
    assert verticall.distance.get_interpolated_median([0.0, 0.0, 1.0, 0.0]) == pytest.approx(2.0)
    assert verticall.distance.get_interpolated_median([0.0, 0.0, 0.0, 1.0]) == pytest.approx(3.0)
    assert verticall.distance.get_interpolated_median([0.5, 0.5, 0.0, 0.0]) == pytest.approx(0.5)
    assert verticall.distance.get_interpolated_median([0.0, 0.5, 0.5, 0.0]) == pytest.approx(1.5)
    assert verticall.distance.get_interpolated_median([0.0, 0.0, 0.5, 0.5]) == pytest.approx(2.5)
    assert 0.000 < verticall.distance.get_interpolated_median([0.7, 0.3, 0.0, 0.0]) < 0.5
    assert 0.125 < verticall.distance.get_interpolated_median([0.3, 0.7, 0.0, 0.0]) < 1.0
    assert 0.250 < verticall.distance.get_interpolated_median([0.0, 0.7, 0.3, 0.0]) < 1.5
    assert 0.375 < verticall.distance.get_interpolated_median([0.0, 0.3, 0.7, 0.0]) < 2.0
    assert 0.500 < verticall.distance.get_interpolated_median([0.0, 0.0, 0.7, 0.3]) < 2.5
    assert 0.625 < verticall.distance.get_interpolated_median([0.0, 0.0, 0.3, 0.7]) < 3.0
    assert (verticall.distance.get_interpolated_median([0.0, 1.000, 0.000, 0.0])
            < verticall.distance.get_interpolated_median([0.0, 0.999, 0.001, 0.0])
            < verticall.distance.get_interpolated_median([0.0, 0.501, 0.499, 0.0])
            < verticall.distance.get_interpolated_median([0.0, 0.500, 0.500, 0.0])
            < verticall.distance.get_interpolated_median([0.0, 0.499, 0.501, 0.0])
            < verticall.distance.get_interpolated_median([0.0, 0.001, 0.990, 0.0])
            < verticall.distance.get_interpolated_median([0.0, 0.000, 1.000, 0.0]))
    assert verticall.distance.get_interpolated_median([0, 0, 1, 0, 10, 9]) == pytest.approx(4.4)
    assert verticall.distance.get_interpolated_median([0, 2, 1, 6, 10, 1]) == pytest.approx(3.6)


def test_get_mode():
    assert verticall.distance.get_mode([1.00, 0.00, 0.00, 0.00]) == 0
    assert verticall.distance.get_mode([0.00, 1.00, 0.00, 0.00]) == 1
    assert verticall.distance.get_mode([0.00, 0.00, 1.00, 0.00]) == 2
    assert verticall.distance.get_mode([0.00, 0.00, 0.00, 1.00]) == 3
    assert verticall.distance.get_mode([0.50, 0.40, 0.00, 0.10]) == 0
    assert verticall.distance.get_mode([0.30, 0.40, 0.10, 0.20]) == 1
    assert verticall.distance.get_mode([0.05, 0.00, 0.90, 0.05]) == 2
    assert verticall.distance.get_mode([0.49, 0.00, 0.00, 0.51]) == 3
    assert verticall.distance.get_mode([0.50, 0.50, 0.00, 0.00]) == pytest.approx(0.5)
    assert verticall.distance.get_mode([0.00, 0.50, 0.50, 0.00]) == pytest.approx(1.5)
    assert verticall.distance.get_mode([0.00, 0.00, 0.50, 0.50]) == pytest.approx(2.5)
    assert verticall.distance.get_mode([0.25, 0.25, 0.25, 0.25]) == pytest.approx(1.5)
    assert verticall.distance.get_mode([0.40, 0.40, 0.00, 0.20]) == pytest.approx(0.5)
    assert verticall.distance.get_mode([0.26, 0.24, 0.26, 0.24]) == pytest.approx(1.0)
    assert verticall.distance.get_mode([0.24, 0.26, 0.24, 0.26]) == pytest.approx(2.0)


def test_interpolate():
    assert verticall.distance.interpolate(0.0, 0.1, 0.0) == pytest.approx(0.0)
    assert verticall.distance.interpolate(0.0, 0.5, 0.0) == pytest.approx(0.0)
    assert verticall.distance.interpolate(0.1, 0.5, 0.1) == pytest.approx(0.0)
    assert verticall.distance.interpolate(0.2, 0.5, 0.2) == pytest.approx(0.0)
    assert verticall.distance.interpolate(0.3, 0.3, 0.3) == pytest.approx(0.0)
    assert verticall.distance.interpolate(0.0, 0.1, 0.1) == pytest.approx(0.5)
    assert verticall.distance.interpolate(0.1, 0.4, 0.4) == pytest.approx(0.5)
    assert verticall.distance.interpolate(0.1, 0.1, 0.0) == pytest.approx(-0.5)
    assert verticall.distance.interpolate(0.4, 0.4, 0.1) == pytest.approx(-0.5)
    assert verticall.distance.interpolate(0.2, 0.4, 0.3) == pytest.approx(0.25)
    assert verticall.distance.interpolate(0.3, 0.4, 0.2) == pytest.approx(-0.25)


def test_climb_to_peak_1():
    # Climb to a single peak at position 0.
    assert verticall.distance.climb_to_peak([0.4, 0.3, 0.2, 0.1, 0.0], 0) == 0
    assert verticall.distance.climb_to_peak([0.4, 0.3, 0.2, 0.1, 0.0], 1) == 0
    assert verticall.distance.climb_to_peak([0.4, 0.3, 0.2, 0.1, 0.0], 2) == 0
    assert verticall.distance.climb_to_peak([0.4, 0.3, 0.2, 0.1, 0.0], 3) == 0
    assert verticall.distance.climb_to_peak([0.4, 0.3, 0.2, 0.1, 0.0], 4) == 0


def test_climb_to_peak_2():
    # Climb to a single peak at position 2.
    assert verticall.distance.climb_to_peak([0.1, 0.3, 0.4, 0.2, 0.0], 0) == 2
    assert verticall.distance.climb_to_peak([0.1, 0.3, 0.4, 0.2, 0.0], 1) == 2
    assert verticall.distance.climb_to_peak([0.1, 0.3, 0.4, 0.2, 0.0], 2) == 2
    assert verticall.distance.climb_to_peak([0.1, 0.3, 0.4, 0.2, 0.0], 3) == 2
    assert verticall.distance.climb_to_peak([0.1, 0.3, 0.4, 0.2, 0.0], 4) == 2


def test_climb_to_peak_3():
    # Two different peaks (positions 1 and 3), so depends on starting position.
    assert verticall.distance.climb_to_peak([0.0, 0.3, 0.2, 0.4, 0.1], 0) == 1
    assert verticall.distance.climb_to_peak([0.0, 0.3, 0.2, 0.4, 0.1], 1) == 1
    assert verticall.distance.climb_to_peak([0.0, 0.3, 0.2, 0.4, 0.1], 2) == 3
    assert verticall.distance.climb_to_peak([0.0, 0.3, 0.2, 0.4, 0.1], 3) == 3
    assert verticall.distance.climb_to_peak([0.0, 0.3, 0.2, 0.4, 0.1], 4) == 3


def test_climb_to_peak_4():
    # When two adjacent positions tie for the peak, the lower one is chosen.
    assert verticall.distance.climb_to_peak([0.35, 0.35, 0.15, 0.10, 0.05], 0) == 0
    assert verticall.distance.climb_to_peak([0.35, 0.35, 0.15, 0.10, 0.05], 1) == 0
    assert verticall.distance.climb_to_peak([0.35, 0.35, 0.15, 0.10, 0.05], 2) == 0
    assert verticall.distance.climb_to_peak([0.35, 0.35, 0.15, 0.10, 0.05], 3) == 0
    assert verticall.distance.climb_to_peak([0.35, 0.35, 0.15, 0.10, 0.05], 4) == 0
    assert verticall.distance.climb_to_peak([0.10, 0.35, 0.35, 0.15, 0.05], 0) == 1
    assert verticall.distance.climb_to_peak([0.10, 0.35, 0.35, 0.15, 0.05], 1) == 1
    assert verticall.distance.climb_to_peak([0.10, 0.35, 0.35, 0.15, 0.05], 2) == 1
    assert verticall.distance.climb_to_peak([0.10, 0.35, 0.35, 0.15, 0.05], 3) == 1
    assert verticall.distance.climb_to_peak([0.10, 0.35, 0.35, 0.15, 0.05], 4) == 1


def test_find_peaks_1():
    # Easy cases with single-point peaks.
    assert verticall.distance.find_peaks([1.0, 0.0, 0.0, 0.0]) == [0]
    assert verticall.distance.find_peaks([0.5, 0.0, 0.5, 0.0]) == [0, 2]
    assert verticall.distance.find_peaks([0.5, 0.0, 0.0, 0.5]) == [0, 3]
    assert verticall.distance.find_peaks([0.5, 0.0, 0.5, 0.0, 0.5]) == [0, 2, 4]
    assert verticall.distance.find_peaks([0.1, 0.2, 0.1, 0.5, 0.1]) == [1, 3]
    assert verticall.distance.find_peaks([0.0, 0.1, 0.1, 0.1, 0.2, 0.2, 0.3, 0.0, 0.0]) == [6]
    assert verticall.distance.find_peaks([0.0, 0.0, 0.3, 0.2, 0.2, 0.1, 0.1, 0.1, 0.0]) == [2]


def test_find_peaks_2():
    # Harder cases with multi-point peaks.
    assert verticall.distance.find_peaks([0.5, 0.5, 0.0, 0.0, 0.0]) == [0]
    assert verticall.distance.find_peaks([0.3, 0.3, 0.3, 0.0, 0.0]) == [1]
    assert verticall.distance.find_peaks([0.0, 0.5, 0.5, 0.0, 0.0]) == [1]
    assert verticall.distance.find_peaks([0.0, 0.3, 0.3, 0.3, 0.0]) == [2]
    assert verticall.distance.find_peaks([0.0, 0.0, 0.0, 0.5, 0.5]) == [3]
    assert verticall.distance.find_peaks([0.0, 0.0, 0.3, 0.3, 0.3]) == [3]
    assert verticall.distance.find_peaks([0.2, 0.2, 0.0, 0.2, 0.2, 0.2, 0.0]) == [0, 4]
    assert verticall.distance.find_peaks([0.0, 0.2, 0.2, 0.0, 0.2, 0.2, 0.2]) == [1, 5]


def test_get_peak_total_mass():
    masses = [0.0, 0.1, 0.2, 0.5, 0.1, 0.1, 0.0]
    assert verticall.distance.get_peak_total_mass(masses, 3) == pytest.approx(1.0)
    masses = [0.1, 0.2, 0.1, 0.0, 0.1, 0.4, 0.1]
    assert verticall.distance.get_peak_total_mass(masses, 1) == pytest.approx(0.4)
    assert verticall.distance.get_peak_total_mass(masses, 5) == pytest.approx(0.6)
    masses = [0.6, 0.1, 0.0, 0.0, 0.0, 0.1, 0.2]
    assert verticall.distance.get_peak_total_mass(masses, 0) == pytest.approx(0.7)
    assert verticall.distance.get_peak_total_mass(masses, 6) == pytest.approx(0.3)


def test_get_sliding_window_count():
    # Test get_window_count() by checking the numbers directly.
    cigars = ['='*1000, '='*100, '='*10]
    assert verticall.distance.get_sliding_window_count(cigars, 1000, 100) == 1
    assert verticall.distance.get_sliding_window_count(cigars, 500, 100) == 6
    assert verticall.distance.get_sliding_window_count(cigars, 100, 100) == 11
    assert verticall.distance.get_sliding_window_count(cigars, 100, 10) == 92
    assert verticall.distance.get_sliding_window_count(cigars, 10, 10) == 111


def test_find_local_minimum_to_right():
    masses = [0.25, 0.20, 0.10, 0.20, 0.25]
    assert verticall.distance.find_local_minimum_to_right(masses, 0) == 2
    masses = [0.25, 0.21, 0.19, 0.10, 0.25]
    assert verticall.distance.find_local_minimum_to_right(masses, 0) == 3
    masses = [0.25, 0.10, 0.19, 0.21, 0.25]
    assert verticall.distance.find_local_minimum_to_right(masses, 0) == 1
    masses = [0.26, 0.24, 0.21, 0.19, 0.10]
    assert verticall.distance.find_local_minimum_to_right(masses, 0) is None
    masses = [0.10, 0.19, 0.21, 0.24, 0.26]
    assert verticall.distance.find_local_minimum_to_right(masses, 4) is None


def test_find_local_minimum_to_left():
    masses = [0.25, 0.20, 0.10, 0.20, 0.25]
    assert verticall.distance.find_local_minimum_to_left(masses, 4) == 2
    masses = [0.25, 0.21, 0.19, 0.10, 0.25]
    assert verticall.distance.find_local_minimum_to_left(masses, 4) == 3
    masses = [0.25, 0.10, 0.19, 0.21, 0.25]
    assert verticall.distance.find_local_minimum_to_left(masses, 4) == 1
    masses = [0.10, 0.19, 0.21, 0.24, 0.26]
    assert verticall.distance.find_local_minimum_to_left(masses, 4) is None
    masses = [0.26, 0.24, 0.21, 0.19, 0.10]
    assert verticall.distance.find_local_minimum_to_left(masses, 0) is None


def test_find_local_maximum_to_right():
    masses = [0.25, 0.10, 0.20, 0.25, 0.20]
    assert verticall.distance.find_local_maximum_to_right(masses, 1) == 3
    masses = [0.25, 0.10, 0.20, 0.25, 0.20]
    assert verticall.distance.find_local_maximum_to_right(masses, 2) == 3
    masses = [0.20, 0.25, 0.10, 0.20, 0.25]
    assert verticall.distance.find_local_maximum_to_right(masses, 2) is None
    assert verticall.distance.find_local_maximum_to_right(masses, 4) is None


def test_find_local_maximum_to_left():
    masses = [0.25, 0.10, 0.20, 0.25, 0.20]
    assert verticall.distance.find_local_maximum_to_left(masses, 1) is None
    assert verticall.distance.find_local_maximum_to_left(masses, 0) is None
    masses = [0.20, 0.25, 0.10, 0.20, 0.25]
    assert verticall.distance.find_local_maximum_to_left(masses, 2) == 1


def test_get_epanechnikov_weight():
    assert verticall.distance.get_epanechnikov_weight(0.0, 0.0) == pytest.approx(1.0)
    assert verticall.distance.get_epanechnikov_weight(0.0, 0.5) == pytest.approx(0.0)
    assert verticall.distance.get_epanechnikov_weight(0.0, -0.5) == pytest.approx(0.0)
    assert verticall.distance.get_epanechnikov_weight(0.0, 5.0) == pytest.approx(0.0)
    assert verticall.distance.get_epanechnikov_weight(0.0, -5.0) == pytest.approx(0.0)

    assert verticall.distance.get_epanechnikov_weight(1.0, 0.0) == pytest.approx(1.0)
    assert verticall.distance.get_epanechnikov_weight(1.0, 0.5) == pytest.approx(0.75)
    assert verticall.distance.get_epanechnikov_weight(1.0, 1.0) == pytest.approx(0.0)
    assert verticall.distance.get_epanechnikov_weight(1.0, 5.0) == pytest.approx(0.0)
    assert verticall.distance.get_epanechnikov_weight(1.0, -0.5) == pytest.approx(0.75)
    assert verticall.distance.get_epanechnikov_weight(1.0, -1.0) == pytest.approx(0.0)
    assert verticall.distance.get_epanechnikov_weight(1.0, -5.0) == pytest.approx(0.0)

    assert verticall.distance.get_epanechnikov_weight(5.0, 0.0) == pytest.approx(1.0)
    assert verticall.distance.get_epanechnikov_weight(5.0, 2.5) == pytest.approx(0.75)
    assert verticall.distance.get_epanechnikov_weight(5.0, 5.0) == pytest.approx(0.0)
    assert verticall.distance.get_epanechnikov_weight(5.0, 10.0) == pytest.approx(0.0)
    assert verticall.distance.get_epanechnikov_weight(5.0, -2.5) == pytest.approx(0.75)
    assert verticall.distance.get_epanechnikov_weight(5.0, -5.0) == pytest.approx(0.0)
    assert verticall.distance.get_epanechnikov_weight(5.0, -10.0) == pytest.approx(0.0)


def test_get_smoothed_mass():
    masses = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    assert verticall.distance.get_smoothed_mass(masses, 0, 0) == pytest.approx(0.1)
    assert verticall.distance.get_smoothed_mass(masses, 0, 1) == pytest.approx(0.1)
    assert verticall.distance.get_smoothed_mass(masses, 0, 5) == pytest.approx(0.1)
    assert verticall.distance.get_smoothed_mass(masses, 3, 0) == pytest.approx(0.1)
    assert verticall.distance.get_smoothed_mass(masses, 3, 1) == pytest.approx(0.1)
    assert verticall.distance.get_smoothed_mass(masses, 3, 5) == pytest.approx(0.1)


def test_get_peak_distance_1():
    masses = [0.0, 0.1, 0.25, 0.1, 0.0, 0.0, 0.1, 0.35, 0.1]
    window_size = 1000
    secondary_ratio = 1.0  # only return the primary result
    mass_peaks, results, _ = \
        verticall.distance.get_peak_distance(masses, window_size, secondary_ratio)
    assert mass_peaks == '0.002000000,0.007000000'
    assert len(results) == 1
    mass, result_level, peak_distance, thresholds = results[0]
    assert result_level == 'primary'
    assert peak_distance == pytest.approx(0.007000000)
    assert mass == pytest.approx(0.55)


def test_get_peak_distance_2():
    masses = [0.0, 0.1, 0.25, 0.1, 0.0, 0.0, 0.1, 0.35, 0.1]
    window_size = 1000
    secondary_ratio = 0.8  # low enough to also return the secondary peak
    mass_peaks, results, _ = \
        verticall.distance.get_peak_distance(masses, window_size, secondary_ratio)
    assert mass_peaks == '0.002000000,0.007000000'
    assert len(results) == 2
    mass, result_level, peak_distance, thresholds = results[0]
    assert result_level == 'primary'
    assert peak_distance == pytest.approx(0.007000000)
    assert mass == pytest.approx(0.55)
    mass, result_level, peak_distance, thresholds = results[1]
    assert result_level == 'secondary'
    assert peak_distance == pytest.approx(0.002000000)
    assert mass == pytest.approx(0.45)


def test_smooth_distribution():
    masses = [0.0000, 0.0006, 0.0009, 0.0012, 0.0017, 0.0023, 0.0033, 0.0044, 0.0056, 0.0059,
              0.0084, 0.0110, 0.0127, 0.0170, 0.0194, 0.0241, 0.0298, 0.0257, 0.0386, 0.0378,
              0.0448, 0.0401, 0.0536, 0.0429, 0.0447, 0.0503, 0.0508, 0.0477, 0.0533, 0.0401,
              0.0481, 0.0360, 0.0341, 0.0280, 0.0224, 0.0199, 0.0186, 0.0166, 0.0145, 0.0092,
              0.0097, 0.0060, 0.0060, 0.0040, 0.0026, 0.0023, 0.0014, 0.0010, 0.0009, 0.0000]
    smoothed_01 = verticall.distance.smooth_distribution(masses, 0.1)
    smoothed_03 = verticall.distance.smooth_distribution(masses, 0.3)
    smoothed_05 = verticall.distance.smooth_distribution(masses, 0.5)

    def delta(data):
        return sum(abs(data[i + 1] - data[i]) for i in range(len(data) - 1))

    assert delta(masses) > delta(smoothed_01) > delta(smoothed_03) > delta(smoothed_05)

    assert verticall.distance.smooth_distribution([], 0.1) == []
    assert verticall.distance.smooth_distribution([], 0.5) == []
