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

import phylo.alignment


def test_get_expanded_cigar():
    assert phylo.alignment.get_expanded_cigar('5=') == '====='
    assert phylo.alignment.get_expanded_cigar('3=1I4=2D2=1X4=') == '===I====DD==X===='
    assert phylo.alignment.get_expanded_cigar('') == ''


def test_overlaps_1():
    # Completely overlapped - cover the same query range.
    a = phylo.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                  'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = phylo.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                  'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 10, 20]:
        assert a.overlaps(b, allowed_overlap)
        assert b.overlaps(a, allowed_overlap)


def test_overlaps_2():
    # Not at all overlapped - far apart on the query.
    a = phylo.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                  'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = phylo.alignment.Alignment('A\t1000\t450\t550\t+\t'
                                  'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 10, 20]:
        assert not a.overlaps(b, allowed_overlap)
        assert not b.overlaps(a, allowed_overlap)


def test_overlaps_3():
    # Adjacent on the query but not overlapped.
    a = phylo.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                  'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = phylo.alignment.Alignment('A\t1000\t150\t250\t+\t'
                                  'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 10, 20]:
        assert not a.overlaps(b, allowed_overlap)
        assert not b.overlaps(a, allowed_overlap)
    for allowed_overlap in [-1, -10, -20]:
        assert a.overlaps(b, allowed_overlap)
        assert b.overlaps(a, allowed_overlap)


def test_overlaps_4():
    # Overlap by 5 bp on the query.
    a = phylo.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                  'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = phylo.alignment.Alignment('A\t1000\t145\t245\t+\t'
                                  'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 1, 2, 3, 4]:
        assert a.overlaps(b, allowed_overlap)
        assert b.overlaps(a, allowed_overlap)
    for allowed_overlap in [5, 10, 100]:
        assert not a.overlaps(b, allowed_overlap)
        assert not b.overlaps(a, allowed_overlap)
