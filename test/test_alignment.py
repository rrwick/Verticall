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

import verticall.alignment
import verticall.paint


def test_get_expanded_cigar():
    assert verticall.alignment.get_expanded_cigar('5=') == '====='
    assert verticall.alignment.get_expanded_cigar('3=1I4=2D2=1X4=') == '===I====DD==X===='
    assert verticall.alignment.get_expanded_cigar('') == ''


def test_overlaps_1():
    # Completely overlapped - cover the same query range.
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 10, 20]:
        assert a.overlaps(b, allowed_overlap)
        assert b.overlaps(a, allowed_overlap)


def test_overlaps_2():
    # Not at all overlapped - far apart on the query.
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('A\t1000\t450\t550\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 10, 20]:
        assert not a.overlaps(b, allowed_overlap)
        assert not b.overlaps(a, allowed_overlap)


def test_overlaps_3():
    # Adjacent on the query but not overlapped.
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('A\t1000\t150\t250\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 10, 20]:
        assert not a.overlaps(b, allowed_overlap)
        assert not b.overlaps(a, allowed_overlap)
    for allowed_overlap in [-1, -10, -20]:
        assert a.overlaps(b, allowed_overlap)
        assert b.overlaps(a, allowed_overlap)


def test_overlaps_4():
    # Overlap by 5 bp on the query.
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('A\t1000\t145\t245\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 1, 2, 3, 4]:
        assert a.overlaps(b, allowed_overlap)
        assert b.overlaps(a, allowed_overlap)
    for allowed_overlap in [5, 10, 100]:
        assert not a.overlaps(b, allowed_overlap)
        assert not b.overlaps(a, allowed_overlap)


def test_swap_insertions_and_deletions():
    assert verticall.alignment.swap_insertions_and_deletions('==========') == '=========='
    assert verticall.alignment.swap_insertions_and_deletions('==I===II==') == '==D===DD=='
    assert verticall.alignment.swap_insertions_and_deletions('=DDD==D=D=') == '=III==I=I='
    assert verticall.alignment.swap_insertions_and_deletions('=IDI=DD=I=') == '=DID=II=D='


def test_remove_indels():
    assert verticall.alignment.remove_indels('========================') == \
           '========================'
    assert verticall.alignment.remove_indels('====X===XXX=====XX=XX===') == \
           '====X===XXX=====XX=XX==='
    assert verticall.alignment.remove_indels('=======IIIII====II===I==') == '================'
    assert verticall.alignment.remove_indels('===DDD===D====DD===DD===') == '================'
    assert verticall.alignment.remove_indels('==X==II===XXXII=IDD==X==') == '==X=====XXX===X=='

    assert verticall.alignment.remove_indels('=====', [4, 5, 6, 7, 8]) == ('=====', [4, 5, 6, 7, 8])
    assert verticall.alignment.remove_indels('==II=', [4, 5, 6, 7, 8]) == ('===', [4, 5, 8])
    assert verticall.alignment.remove_indels('=D=D=', [4, 5, 6, 7, 8]) == ('===', [4, 6, 8])


def test_compress_indels():
    assert verticall.alignment.compress_indels('========================') == \
           '========================'
    assert verticall.alignment.compress_indels('====X===XXX=====XX=XX===') == \
           '====X===XXX=====XX=XX==='
    assert verticall.alignment.compress_indels('=======IIIII====II===I==') == '=======I====I===I=='
    assert verticall.alignment.compress_indels('===DDD===D====DD===DD===') == '===D===D====D===D==='
    assert verticall.alignment.compress_indels('==X==II===XXXII=IDD==X==') == \
           '==X==I===XXXI=ID==X=='

    assert verticall.alignment.compress_indels('=====', [4, 5, 6, 7, 8]) == \
           ('=====', [4, 5, 6, 7, 8])
    assert verticall.alignment.compress_indels('==II=', [4, 5, 6, 7, 8]) == ('==I=', [4, 5, 7, 8])
    assert verticall.alignment.compress_indels('=D=D=', [4, 5, 6, 7, 8]) == \
           ('=D=D=', [4, 5, 6, 7, 8])
    assert verticall.alignment.compress_indels('=DDD=', [4, 5, 6, 7, 8]) == ('=D=', [4, 7, 8])


def test_cigar_to_contig_pos():
    assert verticall.alignment.cigar_to_contig_pos('=====', 0, 5) == [0, 1, 2, 3, 4]
    assert verticall.alignment.cigar_to_contig_pos('=X=X=', 5, 10) == [5, 6, 7, 8, 9]
    assert verticall.alignment.cigar_to_contig_pos('==I==', 2, 7) == [2, 3, 4, 5, 6]
    assert verticall.alignment.cigar_to_contig_pos('==D==', 3, 7) == [3, 4, 5, 5, 6]
    assert verticall.alignment.cigar_to_contig_pos('=D===', 3, 7) == [3, 4, 4, 5, 6]

    assert verticall.alignment.cigar_to_contig_pos('=====', 0, 5, '-') == [4, 3, 2, 1, 0]
    assert verticall.alignment.cigar_to_contig_pos('=X=X=', 5, 10, '-') == [9, 8, 7, 6, 5]
    assert verticall.alignment.cigar_to_contig_pos('==I==', 2, 7, '-') == [6, 5, 4, 3, 2]
    assert verticall.alignment.cigar_to_contig_pos('==D==', 3, 7, '-') == [6, 5, 5, 4, 3]
    assert verticall.alignment.cigar_to_contig_pos('=D===', 3, 7, '-') == [6, 5, 4, 4, 3]


def test_get_difference_count_1():
    assert verticall.alignment.get_difference_count('==================================') == 0
    assert verticall.alignment.get_difference_count('==========X=======================') == 1
    assert verticall.alignment.get_difference_count('====D============D========D=======') == 3
    assert verticall.alignment.get_difference_count('===========I==========I===========') == 2
    assert verticall.alignment.get_difference_count('======D==X===I===XX====I=====D====') == 7


def test_get_difference_count_2():
    cigar = '=====XX===X===DDD===I======IIII===D==X====='
    assert verticall.alignment.get_difference_count(verticall.alignment.compress_indels(cigar)) == 8


def test_find_ambiguous_runs():
    h = verticall.paint.Paint.HORIZONTAL
    v = verticall.paint.Paint.VERTICAL
    a = verticall.paint.Paint.AMBIGUOUS
    assert verticall.alignment.find_ambiguous_runs([h, h, h, h, h, h]) == []
    assert verticall.alignment.find_ambiguous_runs([v, v, v, v, v, v]) == []
    assert verticall.alignment.find_ambiguous_runs([h, h, a, h, h, h]) == [(2, 3)]
    assert verticall.alignment.find_ambiguous_runs([v, v, v, a, a, v]) == [(3, 5)]
    assert verticall.alignment.find_ambiguous_runs([a, a, a, h, a, h]) == [(0, 3), (4, 5)]
    assert verticall.alignment.find_ambiguous_runs([a, v, a, a, a, a]) == [(0, 1), (2, 6)]
    assert verticall.alignment.find_ambiguous_runs([a, a, a, a, a, a]) == [(0, 6)]


def test_remove_ambiguous():
    h = verticall.paint.Paint.HORIZONTAL
    v = verticall.paint.Paint.VERTICAL
    a = verticall.paint.Paint.AMBIGUOUS
    assert verticall.alignment.remove_ambiguous([h, h, h, h, h, h]) == [h, h, h, h, h, h]
    assert verticall.alignment.remove_ambiguous([h, h, a, a, a, h]) == [h, h, h, h, h, h]
    assert verticall.alignment.remove_ambiguous([h, h, a, h, a, h]) == [h, h, h, h, h, h]
    assert verticall.alignment.remove_ambiguous([a, h, h, h, h, a]) == [h, h, h, h, h, h]
    assert verticall.alignment.remove_ambiguous([v, v, v, v, v, v]) == [v, v, v, v, v, v]
    assert verticall.alignment.remove_ambiguous([v, v, a, a, a, v]) == [v, v, v, v, v, v]
    assert verticall.alignment.remove_ambiguous([v, v, a, v, a, v]) == [v, v, v, v, v, v]
    assert verticall.alignment.remove_ambiguous([a, v, v, v, v, a]) == [v, v, v, v, v, v]
    assert verticall.alignment.remove_ambiguous([h, h, h, v, v, v]) == [h, h, h, v, v, v]
    assert verticall.alignment.remove_ambiguous([a, h, h, v, v, a]) == [h, h, h, v, v, v]
    assert verticall.alignment.remove_ambiguous([h, a, h, v, v, a]) == [h, h, h, v, v, v]
    assert verticall.alignment.remove_ambiguous([h, h, a, a, v, v]) == [h, h, h, h, v, v]
    assert verticall.alignment.remove_ambiguous([h, a, a, a, a, v]) == [h, h, h, h, h, v]
    assert verticall.alignment.remove_ambiguous([a, a, a, a, a, a]) == [h, h, h, h, h, h]
