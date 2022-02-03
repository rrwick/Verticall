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

import phylo.align


def test_compress_indels():
    assert phylo.align.compress_indels('========================') == '========================'
    assert phylo.align.compress_indels('====X===XXX=====XX=XX===') == '====X===XXX=====XX=XX==='
    assert phylo.align.compress_indels('=======IIIII====II===I==') == '=======I====I===I=='
    assert phylo.align.compress_indels('===DDD===D====DD===DD===') == '===D===D====D===D==='
    assert phylo.align.compress_indels('==X==II===XXXII=IDD==X==') == '==X==I===XXXI=ID==X=='


def test_get_difference_count_1():
    assert phylo.align.get_difference_count('==================================') == 0
    assert phylo.align.get_difference_count('==========X=======================') == 1
    assert phylo.align.get_difference_count('====D============D========D=======') == 3
    assert phylo.align.get_difference_count('===========I==========I===========') == 2
    assert phylo.align.get_difference_count('======D==X===I===XX====I=====D====') == 7


def test_get_difference_count_2():
    cigar = '=====XX===X===DDD===I======IIII===D==X====='
    assert phylo.align.get_difference_count(phylo.align.compress_indels(cigar)) == 8
