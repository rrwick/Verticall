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


def test_compress_indels():
    assert phylo.alignment.compress_indels('========================') == '========================'
    assert phylo.alignment.compress_indels('====X===XXX=====XX=XX===') == '====X===XXX=====XX=XX==='
    assert phylo.alignment.compress_indels('=======IIIII====II===I==') == '=======I====I===I=='
    assert phylo.alignment.compress_indels('===DDD===D====DD===DD===') == '===D===D====D===D==='
    assert phylo.alignment.compress_indels('==X==II===XXXII=IDD==X==') == '==X==I===XXXI=ID==X=='


def test_remove_indels():
    assert phylo.alignment.remove_indels('========================') == '========================'
    assert phylo.alignment.remove_indels('====X===XXX=====XX=XX===') == '====X===XXX=====XX=XX==='
    assert phylo.alignment.remove_indels('=======IIIII====II===I==') == '================'
    assert phylo.alignment.remove_indels('===DDD===D====DD===DD===') == '================'
    assert phylo.alignment.remove_indels('==X==II===XXXII=IDD==X==') == '==X=====XXX===X=='


def test_get_difference_count_1():
    assert phylo.alignment.get_difference_count('==================================') == 0
    assert phylo.alignment.get_difference_count('==========X=======================') == 1
    assert phylo.alignment.get_difference_count('====D============D========D=======') == 3
    assert phylo.alignment.get_difference_count('===========I==========I===========') == 2
    assert phylo.alignment.get_difference_count('======D==X===I===XX====I=====D====') == 7


def test_get_difference_count_2():
    cigar = '=====XX===X===DDD===I======IIII===D==X====='
    assert phylo.alignment.get_difference_count(phylo.alignment.compress_indels(cigar)) == 8


def test_get_window_count_1():
    # Test get_window_count() by checking the numbers directly.
    cigars = ['='*1000, '='*100, '='*10]
    assert phylo.alignment.get_window_count(cigars, 1000, 100) == 1
    assert phylo.alignment.get_window_count(cigars, 500, 100) == 6
    assert phylo.alignment.get_window_count(cigars, 100, 100) == 11
    assert phylo.alignment.get_window_count(cigars, 100, 10) == 92
    assert phylo.alignment.get_window_count(cigars, 10, 10) == 111


def test_get_window_count_2():
    # Test get_window_count() by comparing the numbers to the results of get_distances().
    for cigars in [['='*56743, '='*4324, '='*19], ['='*4565, '='*343, '='*12], ['='*9867, '='*789]]:
        for window_size, window_step in [(1000, 100), (456, 12), (34234, 232), (2938, 23)]:
            assert phylo.alignment.get_window_count(cigars, window_size, window_step) == \
                   len(phylo.alignment.get_distances(cigars, window_size, window_step)[0])
