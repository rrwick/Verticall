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


def test_swap_insertions_and_deletions():
    assert phylo.alignment.swap_insertions_and_deletions('==========') == '=========='
    assert phylo.alignment.swap_insertions_and_deletions('==I===II==') == '==D===DD=='
    assert phylo.alignment.swap_insertions_and_deletions('=DDD==D=D=') == '=III==I=I='
    assert phylo.alignment.swap_insertions_and_deletions('=IDI=DD=I=') == '=DID=II=D='


def test_remove_indels():
    assert phylo.alignment.remove_indels('========================') == '========================'
    assert phylo.alignment.remove_indels('====X===XXX=====XX=XX===') == '====X===XXX=====XX=XX==='
    assert phylo.alignment.remove_indels('=======IIIII====II===I==') == '================'
    assert phylo.alignment.remove_indels('===DDD===D====DD===DD===') == '================'
    assert phylo.alignment.remove_indels('==X==II===XXXII=IDD==X==') == '==X=====XXX===X=='

    assert phylo.alignment.remove_indels('=====', [4, 5, 6, 7, 8]) == ('=====', [4, 5, 6, 7, 8])
    assert phylo.alignment.remove_indels('==II=', [4, 5, 6, 7, 8]) == ('===', [4, 5, 8])
    assert phylo.alignment.remove_indels('=D=D=', [4, 5, 6, 7, 8]) == ('===', [4, 6, 8])


def test_compress_indels():
    assert phylo.alignment.compress_indels('========================') == '========================'
    assert phylo.alignment.compress_indels('====X===XXX=====XX=XX===') == '====X===XXX=====XX=XX==='
    assert phylo.alignment.compress_indels('=======IIIII====II===I==') == '=======I====I===I=='
    assert phylo.alignment.compress_indels('===DDD===D====DD===DD===') == '===D===D====D===D==='
    assert phylo.alignment.compress_indels('==X==II===XXXII=IDD==X==') == '==X==I===XXXI=ID==X=='

    assert phylo.alignment.compress_indels('=====', [4, 5, 6, 7, 8]) == ('=====', [4, 5, 6, 7, 8])
    assert phylo.alignment.compress_indels('==II=', [4, 5, 6, 7, 8]) == ('==I=', [4, 5, 7, 8])
    assert phylo.alignment.compress_indels('=D=D=', [4, 5, 6, 7, 8]) == ('=D=D=', [4, 5, 6, 7, 8])
    assert phylo.alignment.compress_indels('=DDD=', [4, 5, 6, 7, 8]) == ('=D=', [4, 7, 8])


def test_cigar_to_contig_positions():
    assert phylo.alignment.cigar_to_contig_positions('=====', 0, 5) == [0, 1, 2, 3, 4]
    assert phylo.alignment.cigar_to_contig_positions('=X=X=', 5, 10) == [5, 6, 7, 8, 9]
    assert phylo.alignment.cigar_to_contig_positions('==I==', 2, 7) == [2, 3, 4, 5, 6]
    assert phylo.alignment.cigar_to_contig_positions('==D==', 3, 7) == [3, 4, 5, 5, 6]
    assert phylo.alignment.cigar_to_contig_positions('=D===', 3, 7) == [3, 4, 4, 5, 6]
