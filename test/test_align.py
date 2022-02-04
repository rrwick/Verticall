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

import phylo.align


def test_compress_indels():
    assert phylo.align.compress_indels('========================') == '========================'
    assert phylo.align.compress_indels('====X===XXX=====XX=XX===') == '====X===XXX=====XX=XX==='
    assert phylo.align.compress_indels('=======IIIII====II===I==') == '=======I====I===I=='
    assert phylo.align.compress_indels('===DDD===D====DD===DD===') == '===D===D====D===D==='
    assert phylo.align.compress_indels('==X==II===XXXII=IDD==X==') == '==X==I===XXXI=ID==X=='


def test_remove_indels():
    assert phylo.align.remove_indels('========================') == '========================'
    assert phylo.align.remove_indels('====X===XXX=====XX=XX===') == '====X===XXX=====XX=XX==='
    assert phylo.align.remove_indels('=======IIIII====II===I==') == '================'
    assert phylo.align.remove_indels('===DDD===D====DD===DD===') == '================'
    assert phylo.align.remove_indels('==X==II===XXXII=IDD==X==') == '==X=====XXX===X=='


def test_get_difference_count_1():
    assert phylo.align.get_difference_count('==================================') == 0
    assert phylo.align.get_difference_count('==========X=======================') == 1
    assert phylo.align.get_difference_count('====D============D========D=======') == 3
    assert phylo.align.get_difference_count('===========I==========I===========') == 2
    assert phylo.align.get_difference_count('======D==X===I===XX====I=====D====') == 7


def test_get_difference_count_2():
    cigar = '=====XX===X===DDD===I======IIII===D==X====='
    assert phylo.align.get_difference_count(phylo.align.compress_indels(cigar)) == 8


def test_get_window_count_1():
    # Test get_window_count() by checking the numbers directly.
    cigars = ['='*1000, '='*100, '='*10]
    assert phylo.align.get_window_count(cigars, 1000, 100) == 1
    assert phylo.align.get_window_count(cigars, 500, 100) == 6
    assert phylo.align.get_window_count(cigars, 100, 100) == 11
    assert phylo.align.get_window_count(cigars, 100, 10) == 92
    assert phylo.align.get_window_count(cigars, 10, 10) == 111


def test_get_window_count_2():
    # Test get_window_count() by comparing the numbers to the results of get_distances().
    for cigars in [['='*56743, '='*4324, '='*19], ['='*4565, '='*343, '='*12], ['='*9867, '='*789]]:
        for window_size, window_step in [(1000, 100), (456, 12), (34234, 232), (2938, 23)]:
            assert phylo.align.get_window_count(cigars, window_size, window_step) == \
                   len(phylo.align.get_distances(cigars, window_size, window_step)[0])
