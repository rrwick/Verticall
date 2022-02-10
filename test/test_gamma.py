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
from scipy.stats import gamma

import phylo.gamma


def test_get_mean():
    assert phylo.gamma.get_mean([1.00, 0.00, 0.00, 0.00]) == pytest.approx(0.0)
    assert phylo.gamma.get_mean([0.00, 1.00, 0.00, 0.00]) == pytest.approx(1.0)
    assert phylo.gamma.get_mean([0.00, 0.00, 1.00, 0.00]) == pytest.approx(2.0)
    assert phylo.gamma.get_mean([0.00, 0.00, 0.00, 1.00]) == pytest.approx(3.0)
    assert phylo.gamma.get_mean([0.50, 0.50, 0.00, 0.00]) == pytest.approx(0.5)
    assert phylo.gamma.get_mean([0.00, 0.50, 0.50, 0.00]) == pytest.approx(1.5)
    assert phylo.gamma.get_mean([0.00, 0.00, 0.50, 0.50]) == pytest.approx(2.5)
    assert phylo.gamma.get_mean([0.25, 0.25, 0.25, 0.25]) == pytest.approx(1.5)
    assert phylo.gamma.get_mean([0.10, 0.20, 0.30, 0.40]) == pytest.approx(2.0)
    assert phylo.gamma.get_mean([0.40, 0.30, 0.20, 0.10]) == pytest.approx(1.0)


def test_get_variance():
    assert phylo.gamma.get_variance([1.0, 0.0, 0.0, 0.0]) == pytest.approx(0.0)
    assert phylo.gamma.get_variance([0.0, 1.0, 0.0, 0.0]) == pytest.approx(0.0)
    assert phylo.gamma.get_variance([0.0, 0.0, 1.0, 0.0]) == pytest.approx(0.0)
    assert phylo.gamma.get_variance([0.0, 0.0, 0.0, 1.0]) == pytest.approx(0.0)
    assert phylo.gamma.get_variance([0.2, 0.2, 0.2, 0.2, 0.2]) == pytest.approx(2.0)
    assert phylo.gamma.get_variance([0.1, 0.2, 0.4, 0.3]) == pytest.approx(0.89)
    assert phylo.gamma.get_variance([0.16, 0.53, 0.2, 0.08, 0.03]) == pytest.approx(0.8659)
    assert phylo.gamma.get_variance([0, 1, 1, 1, 1, 1, 1]) == pytest.approx(35.0 / 12.0)


def test_get_shape_and_scale():
    for mean, variance in [(5.0, 5.0), (1.0, 5.0), (5.0, 1.0), (5.0, 10.0), (10.0, 5.0)]:
        shape, scale = phylo.gamma.get_shape_and_scale(mean, variance)
        assert gamma.mean(shape, scale=scale) == pytest.approx(mean)
        assert gamma.var(shape, scale=scale) == pytest.approx(variance)
