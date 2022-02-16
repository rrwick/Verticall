"""
This module contains some tests for XXXXXXXXX. To run them, execute `pytest` from the root
XXXXXXXXX directory.

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

import pytest

import phylo.distance


def test_get_mean():
    assert phylo.distance.get_mean([1.00, 0.00, 0.00, 0.00]) == pytest.approx(0.0)
    assert phylo.distance.get_mean([0.00, 1.00, 0.00, 0.00]) == pytest.approx(1.0)
    assert phylo.distance.get_mean([0.00, 0.00, 1.00, 0.00]) == pytest.approx(2.0)
    assert phylo.distance.get_mean([0.00, 0.00, 0.00, 1.00]) == pytest.approx(3.0)
    assert phylo.distance.get_mean([0.50, 0.50, 0.00, 0.00]) == pytest.approx(0.5)
    assert phylo.distance.get_mean([0.00, 0.50, 0.50, 0.00]) == pytest.approx(1.5)
    assert phylo.distance.get_mean([0.00, 0.00, 0.50, 0.50]) == pytest.approx(2.5)
    assert phylo.distance.get_mean([0.25, 0.25, 0.25, 0.25]) == pytest.approx(1.5)
    assert phylo.distance.get_mean([0.10, 0.20, 0.30, 0.40]) == pytest.approx(2.0)
    assert phylo.distance.get_mean([0.40, 0.30, 0.20, 0.10]) == pytest.approx(1.0)


def test_get_median():
    assert phylo.distance.get_median([1.0, 0.0, 0.0, 0.0]) == 0
    assert phylo.distance.get_median([0.0, 1.0, 0.0, 0.0]) == 1
    assert phylo.distance.get_median([0.0, 0.0, 1.0, 0.0]) == 2
    assert phylo.distance.get_median([0.0, 0.0, 0.0, 1.0]) == 3
    assert phylo.distance.get_median([0.6, 0.0, 0.0, 0.4]) == 0
    assert phylo.distance.get_median([0.4, 0.0, 0.0, 0.6]) == 3
    assert phylo.distance.get_median([0.1, 0.2, 0.3, 0.4]) == 2
    assert phylo.distance.get_median([0.4, 0.3, 0.2, 0.1]) == 1


def test_get_interpolated_median():
    assert phylo.distance.get_interpolated_median([1.0, 0.0, 0.0, 0.0]) == pytest.approx(0.0)
    assert phylo.distance.get_interpolated_median([0.0, 1.0, 0.0, 0.0]) == pytest.approx(1.0)
    assert phylo.distance.get_interpolated_median([0.0, 0.0, 1.0, 0.0]) == pytest.approx(2.0)
    assert phylo.distance.get_interpolated_median([0.0, 0.0, 0.0, 1.0]) == pytest.approx(3.0)
    assert phylo.distance.get_interpolated_median([0.5, 0.5, 0.0, 0.0]) == pytest.approx(0.5)
    assert phylo.distance.get_interpolated_median([0.0, 0.5, 0.5, 0.0]) == pytest.approx(1.5)
    assert phylo.distance.get_interpolated_median([0.0, 0.0, 0.5, 0.5]) == pytest.approx(2.5)
    assert 0.000 < phylo.distance.get_interpolated_median([0.7, 0.3, 0.0, 0.0]) < 0.5
    assert 0.125 < phylo.distance.get_interpolated_median([0.3, 0.7, 0.0, 0.0]) < 1.0
    assert 0.250 < phylo.distance.get_interpolated_median([0.0, 0.7, 0.3, 0.0]) < 1.5
    assert 0.375 < phylo.distance.get_interpolated_median([0.0, 0.3, 0.7, 0.0]) < 2.0
    assert 0.500 < phylo.distance.get_interpolated_median([0.0, 0.0, 0.7, 0.3]) < 2.5
    assert 0.625 < phylo.distance.get_interpolated_median([0.0, 0.0, 0.3, 0.7]) < 3.0
    assert (phylo.distance.get_interpolated_median([0.0, 1.000, 0.000, 0.0])
            < phylo.distance.get_interpolated_median([0.0, 0.999, 0.001, 0.0])
            < phylo.distance.get_interpolated_median([0.0, 0.501, 0.499, 0.0])
            < phylo.distance.get_interpolated_median([0.0, 0.500, 0.500, 0.0])
            < phylo.distance.get_interpolated_median([0.0, 0.499, 0.501, 0.0])
            < phylo.distance.get_interpolated_median([0.0, 0.001, 0.990, 0.0])
            < phylo.distance.get_interpolated_median([0.0, 0.000, 1.000, 0.0]))
    assert phylo.distance.get_interpolated_median([0, 0, 1, 0, 10, 9]) == pytest.approx(4.4)
    assert phylo.distance.get_interpolated_median([0, 2, 1, 6, 10, 1]) == pytest.approx(3.6)


def test_get_mode():
    assert phylo.distance.get_mode([1.00, 0.00, 0.00, 0.00]) == 0
    assert phylo.distance.get_mode([0.00, 1.00, 0.00, 0.00]) == 1
    assert phylo.distance.get_mode([0.00, 0.00, 1.00, 0.00]) == 2
    assert phylo.distance.get_mode([0.00, 0.00, 0.00, 1.00]) == 3
    assert phylo.distance.get_mode([0.50, 0.40, 0.00, 0.10]) == 0
    assert phylo.distance.get_mode([0.30, 0.40, 0.10, 0.20]) == 1
    assert phylo.distance.get_mode([0.05, 0.00, 0.90, 0.05]) == 2
    assert phylo.distance.get_mode([0.49, 0.00, 0.00, 0.51]) == 3
    assert phylo.distance.get_mode([0.50, 0.50, 0.00, 0.00]) == pytest.approx(0.5)
    assert phylo.distance.get_mode([0.00, 0.50, 0.50, 0.00]) == pytest.approx(1.5)
    assert phylo.distance.get_mode([0.00, 0.00, 0.50, 0.50]) == pytest.approx(2.5)
    assert phylo.distance.get_mode([0.25, 0.25, 0.25, 0.25]) == pytest.approx(1.5)
    assert phylo.distance.get_mode([0.40, 0.40, 0.00, 0.20]) == pytest.approx(0.5)
    assert phylo.distance.get_mode([0.26, 0.24, 0.26, 0.24]) == pytest.approx(1.0)
    assert phylo.distance.get_mode([0.24, 0.26, 0.24, 0.26]) == pytest.approx(2.0)


def test_get_top_half_1():
    # Simple cases with all mass in one point.
    assert phylo.distance.get_top_fraction([1.00, 0.00, 0.00, 0.00], 0.5) == (0, 1)
    assert phylo.distance.get_top_fraction([0.00, 1.00, 0.00, 0.00], 0.5) == (1, 2)
    assert phylo.distance.get_top_fraction([0.00, 0.00, 1.00, 0.00], 0.5) == (2, 3)
    assert phylo.distance.get_top_fraction([0.00, 0.00, 0.00, 1.00], 0.5) == (3, 4)


def test_get_top_half_2():
    # Ranges which start on zero will be extended to the right as long as they
    assert phylo.distance.get_top_fraction([0.60, 0.00, 0.00, 0.00, 0.00], 0.5) == (0, 1)
    assert phylo.distance.get_top_fraction([0.60, 0.10, 0.00, 0.00, 0.00], 0.5) == (0, 2)
    assert phylo.distance.get_top_fraction([0.60, 0.10, 0.05, 0.00, 0.00], 0.5) == (0, 3)
    assert phylo.distance.get_top_fraction([0.60, 0.10, 0.10, 0.00, 0.00], 0.5) == (0, 2)
    assert phylo.distance.get_top_fraction([0.60, 0.10, 0.15, 0.00, 0.00], 0.5) == (0, 2)
    assert phylo.distance.get_top_fraction([0.60, 0.10, 0.09, 0.08, 0.07], 0.5) == (0, 5)


def test_get_top_half_3():
    # A tie will expand in both directions.
    assert phylo.distance.get_top_fraction([0.30, 0.40, 0.30, 0.00, 0.00], 0.5) == (0, 3)
    assert phylo.distance.get_top_fraction([0.00, 0.30, 0.40, 0.30, 0.00], 0.5) == (1, 4)
    assert phylo.distance.get_top_fraction([0.00, 0.00, 0.30, 0.40, 0.30], 0.5) == (2, 5)
    assert phylo.distance.get_top_fraction([0.1, 0.15, 0.2, 0.2, 0.15, 0.1, 0.1], 0.5) == (1, 5)
    assert phylo.distance.get_top_fraction([0.1, 0.1, 0.15, 0.2, 0.2, 0.15, 0.1], 0.5) == (2, 6)
    assert phylo.distance.get_top_fraction([0.1, 0.1, 0.1, 0.15, 0.2, 0.2, 0.15], 0.5) == (3, 7)


def test_get_top_half_4():
    # Testing right-skewed and left-skewed distributions
    left_skewed = [0.00, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.06, 0.04, 0.00, 0.00, 0.00]
    right_skewed = [0.00, 0.00, 0.00, 0.04, 0.06, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.00]
    assert phylo.distance.get_top_fraction(left_skewed, 0.5) == (1, 4)
    assert phylo.distance.get_top_fraction(right_skewed, 0.5) == (8, 11)


def test_correct_distances_1():
    sample_names = ['a', 'b']
    distances = {('a', 'a'): 0.0, ('a', 'b'): 0.2,
                 ('b', 'a'): 0.1, ('b', 'b'): 0.0}
    aligned_fractions = {('a', 'a'): 1.0, ('a', 'b'): 0.8,
                         ('b', 'a'): 0.9, ('b', 'b'): 1.0}
    phylo.distance.correct_distances(distances, aligned_fractions, sample_names, {'none'})
    assert distances[('a', 'a')] == pytest.approx(0.0)
    assert distances[('a', 'b')] == pytest.approx(0.2)
    assert distances[('b', 'a')] == pytest.approx(0.1)
    assert distances[('b', 'b')] == pytest.approx(0.0)


def test_correct_distances_2():
    sample_names = ['a', 'b']
    distances = {('a', 'a'): 0.0, ('a', 'b'): 0.2,
                 ('b', 'a'): 0.1, ('b', 'b'): 0.0}
    aligned_fractions = {('a', 'a'): 1.0, ('a', 'b'): 0.8,
                         ('b', 'a'): 0.9, ('b', 'b'): 1.0}
    phylo.distance.correct_distances(distances, aligned_fractions, sample_names, {'jukescantor'})
    assert distances[('a', 'a')] == pytest.approx(0.0)
    assert distances[('a', 'b')] == pytest.approx(0.23261619622788)
    assert distances[('b', 'a')] == pytest.approx(0.107325632730505)
    assert distances[('b', 'b')] == pytest.approx(0.0)


def test_correct_distances_3():
    sample_names = ['a', 'b']
    distances = {('a', 'a'): 0.0, ('a', 'b'): 0.2,
                 ('b', 'a'): 0.1, ('b', 'b'): 0.0}
    aligned_fractions = {('a', 'a'): 1.0, ('a', 'b'): 0.8,
                         ('b', 'a'): 0.9, ('b', 'b'): 1.0}
    phylo.distance.correct_distances(distances, aligned_fractions, sample_names, {'alignedfrac'})
    assert distances[('a', 'a')] == pytest.approx(0.0)
    assert distances[('a', 'b')] == pytest.approx(0.2 / 0.8)
    assert distances[('b', 'a')] == pytest.approx(0.1 / 0.9)
    assert distances[('b', 'b')] == pytest.approx(0.0)


def test_correct_distances_4():
    sample_names = ['a', 'b']
    distances = {('a', 'a'): 0.0, ('a', 'b'): 0.2,
                 ('b', 'a'): 0.1, ('b', 'b'): 0.0}
    aligned_fractions = {('a', 'a'): 1.0, ('a', 'b'): 0.8,
                         ('b', 'a'): 0.9, ('b', 'b'): 1.0}
    phylo.distance.correct_distances(distances, aligned_fractions, sample_names,
                                     {'jukescantor', 'alignedfrac'})
    assert distances[('a', 'a')] == pytest.approx(0.0)
    assert distances[('a', 'b')] == pytest.approx(0.23261619622788 / 0.8)
    assert distances[('b', 'a')] == pytest.approx(0.107325632730505 / 0.9)
    assert distances[('b', 'b')] == pytest.approx(0.0)


def test_make_symmetrical():
    sample_names = ['a', 'b']
    distances = {('a', 'a'): 0.0, ('a', 'b'): 0.2,
                 ('b', 'a'): 0.1, ('b', 'b'): 0.0}
    phylo.distance.make_symmetrical(distances, sample_names)
    assert distances[('a', 'a')] == pytest.approx(0.0)
    assert distances[('a', 'b')] == pytest.approx(0.15)
    assert distances[('b', 'a')] == pytest.approx(0.15)
    assert distances[('b', 'b')] == pytest.approx(0.0)
