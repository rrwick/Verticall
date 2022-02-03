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
If not, see <http://www.gnu.org/licenses/>.
"""

import phylo.intrange


def test_range_1():
    r = phylo.intrange.IntRange()
    assert r.total_length() == 0
    assert str(r) == '[]'


def test_range_2():
    r = phylo.intrange.IntRange()
    r.add_range(0, 10)
    assert r.total_length() == 10
    assert str(r) == '[(0, 10)]'


def test_range_3():
    r = phylo.intrange.IntRange()
    r.add_range(10, 0)
    assert r.total_length() == 10
    assert str(r) == '[(0, 10)]'


def test_range_4():
    r = phylo.intrange.IntRange()
    r.add_range(0, 10)
    r.add_range(20, 30)
    assert r.total_length() == 20
    assert str(r) == '[(0, 10), (20, 30)]'


def test_range_5():
    r = phylo.intrange.IntRange()
    r.add_range(0, 10)
    r.add_range(0, 10)
    r.add_range(0, 10)
    r.add_range(0, 10)
    r.add_range(0, 10)
    assert r.total_length() == 10
    assert str(r) == '[(0, 10)]'


def test_range_6():
    r = phylo.intrange.IntRange()
    r.add_range(0, 10)
    r.add_range(5, 15)
    assert r.total_length() == 15
    assert str(r) == '[(0, 15)]'


def test_overlaps_1():
    r1 = phylo.intrange.IntRange([(0, 10)])
    r2 = phylo.intrange.IntRange([(20, 30)])
    assert not r1.overlaps(r2)
    assert not r2.overlaps(r1)


def test_overlaps_2():
    r1 = phylo.intrange.IntRange([(0, 10)])
    r2 = phylo.intrange.IntRange([(5, 15)])
    assert r1.overlaps(r2)
    assert r2.overlaps(r1)
