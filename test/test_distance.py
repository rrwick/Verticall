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


def test_get_mean_distance():
    assert phylo.distance.get_mean_distance([1.00, 0.00, 0.00, 0.00], 4) == pytest.approx(0.000)
    assert phylo.distance.get_mean_distance([0.00, 1.00, 0.00, 0.00], 4) == pytest.approx(0.250)
    assert phylo.distance.get_mean_distance([0.00, 0.00, 1.00, 0.00], 4) == pytest.approx(0.500)
    assert phylo.distance.get_mean_distance([0.00, 0.00, 0.00, 1.00], 4) == pytest.approx(0.750)
    assert phylo.distance.get_mean_distance([0.50, 0.50, 0.00, 0.00], 4) == pytest.approx(0.125)
    assert phylo.distance.get_mean_distance([0.00, 0.50, 0.50, 0.00], 4) == pytest.approx(0.375)
    assert phylo.distance.get_mean_distance([0.00, 0.00, 0.50, 0.50], 4) == pytest.approx(0.625)
    assert phylo.distance.get_mean_distance([0.25, 0.25, 0.25, 0.25], 4) == pytest.approx(0.375)
    assert phylo.distance.get_mean_distance([0.10, 0.20, 0.30, 0.40], 4) == pytest.approx(0.500)
    assert phylo.distance.get_mean_distance([0.40, 0.30, 0.20, 0.10], 4) == pytest.approx(0.250)


def test_get_median_distance():
    assert phylo.distance.get_median_distance([1.0, 0.0, 0.0, 0.0], 4) == pytest.approx(0.00)
    assert phylo.distance.get_median_distance([0.0, 1.0, 0.0, 0.0], 4) == pytest.approx(0.25)
    assert phylo.distance.get_median_distance([0.0, 0.0, 1.0, 0.0], 4) == pytest.approx(0.50)
    assert phylo.distance.get_median_distance([0.0, 0.0, 0.0, 1.0], 4) == pytest.approx(0.75)
    assert phylo.distance.get_median_distance([0.6, 0.0, 0.0, 0.4], 4) == pytest.approx(0.00)
    assert phylo.distance.get_median_distance([0.4, 0.0, 0.0, 0.6], 4) == pytest.approx(0.75)
    assert phylo.distance.get_median_distance([0.1, 0.2, 0.3, 0.4], 4) == pytest.approx(0.50)
    assert phylo.distance.get_median_distance([0.4, 0.3, 0.2, 0.1], 4) == pytest.approx(0.25)
