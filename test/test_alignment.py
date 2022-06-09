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

import pathlib
import pytest
import tempfile

import verticall.alignment


def test_index_exists():
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir = pathlib.Path(temp_dir)

        assert not verticall.alignment.index_exists(temp_dir, 'absent', True)

        empty_index = temp_dir / 'empty.mmi'
        open(empty_index, 'a').close()
        assert not verticall.alignment.index_exists(temp_dir, 'empty', True)

        good_index = temp_dir / 'good.mmi'
        with open(good_index, 'wt') as f:
            f.write('stuff')
        assert verticall.alignment.index_exists(temp_dir, 'good', True)


def test_get_expanded_cigar():
    assert verticall.alignment.get_expanded_cigar('5=') == '====='
    assert verticall.alignment.get_expanded_cigar('3=1I4=2D2=1X4=') == '===I====DD==X===='
    assert verticall.alignment.get_expanded_cigar('') == ''


def test_bad_paf():
    with pytest.raises(SystemExit) as e:
        verticall.alignment.Alignment('not_a_paf_line')
    assert 'PAF format' in str(e.value)


def test_repr():
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t60\t160\t100\t100\tAS:i:100\tcg:Z:100=')
    assert str(a) == 'A:50-150(+), C:60-160 (100.000%)'


def test_query_covered_bases():
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    assert a.query_covered_bases() == 100


def test_overlaps_on_query_1():
    # Completely overlapped - cover the same query range.
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 10, 20]:
        assert a.overlaps_on_query(b, allowed_overlap)
        assert b.overlaps_on_query(a, allowed_overlap)


def test_overlaps_on_query_2():
    # Not at all overlapped - far apart on the query.
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('A\t1000\t450\t550\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 10, 20]:
        assert not a.overlaps_on_query(b, allowed_overlap)
        assert not b.overlaps_on_query(a, allowed_overlap)


def test_overlaps_on_query_3():
    # Adjacent on the query but not overlapped.
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('A\t1000\t150\t250\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 10, 20]:
        assert not a.overlaps_on_query(b, allowed_overlap)
        assert not b.overlaps_on_query(a, allowed_overlap)
    for allowed_overlap in [-1, -10, -20]:
        assert a.overlaps_on_query(b, allowed_overlap)
        assert b.overlaps_on_query(a, allowed_overlap)


def test_overlaps_on_query_4():
    # Overlap by 5 bp on the query.
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('A\t1000\t145\t245\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 1, 2, 3, 4]:
        assert a.overlaps_on_query(b, allowed_overlap)
        assert b.overlaps_on_query(a, allowed_overlap)
    for allowed_overlap in [5, 10, 100]:
        assert not a.overlaps_on_query(b, allowed_overlap)
        assert not b.overlaps_on_query(a, allowed_overlap)


def test_overlaps_on_query_5():
    # Different query sequence, so no overlap.
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('B\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 10, 20]:
        assert not a.overlaps_on_query(b, allowed_overlap)
        assert not b.overlaps_on_query(a, allowed_overlap)


def test_overlaps_on_target_1():
    # Completely overlapped - cover the same target range.
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 10, 20]:
        assert a.overlaps_on_target(b, allowed_overlap)
        assert b.overlaps_on_target(a, allowed_overlap)


def test_overlaps_on_target_2():
    # Not at all overlapped - far apart on the query.
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t450\t550\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 10, 20]:
        assert not a.overlaps_on_target(b, allowed_overlap)
        assert not b.overlaps_on_target(a, allowed_overlap)


def test_overlaps_on_target_3():
    # Adjacent on the query but not overlapped.
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t150\t250\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 10, 20]:
        assert not a.overlaps_on_target(b, allowed_overlap)
        assert not b.overlaps_on_target(a, allowed_overlap)
    for allowed_overlap in [-1, -10, -20]:
        assert a.overlaps_on_target(b, allowed_overlap)
        assert b.overlaps_on_target(a, allowed_overlap)


def test_overlaps_on_target_4():
    # Overlap by 5 bp on the query.
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t145\t245\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 1, 2, 3, 4]:
        assert a.overlaps_on_target(b, allowed_overlap)
        assert b.overlaps_on_target(a, allowed_overlap)
    for allowed_overlap in [5, 10, 100]:
        assert not a.overlaps_on_target(b, allowed_overlap)
        assert not b.overlaps_on_target(a, allowed_overlap)


def test_overlaps_on_target_5():
    # Different target sequence, so no overlap.
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'B\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 10, 20]:
        assert not a.overlaps_on_target(b, allowed_overlap)
        assert not b.overlaps_on_target(a, allowed_overlap)


def test_overlaps_1():
    # Completely overlapped - cover the same query and target range.
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 10, 20]:
        assert a.overlaps(b, allowed_overlap)
        assert b.overlaps(a, allowed_overlap)


def test_overlaps_2():
    # Completely overlapped in the query (but different targets).
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'B\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 10, 20]:
        assert a.overlaps(b, allowed_overlap)
        assert b.overlaps(a, allowed_overlap)


def test_overlaps_3():
    # Completely overlapped in the targets (but different queries).
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('B\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    for allowed_overlap in [0, 10, 20]:
        assert a.overlaps(b, allowed_overlap)
        assert b.overlaps(a, allowed_overlap)


def test_cull_redundant_alignments_1():
    # Totally different queries and targets, so no culling.
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'B\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('C\t1000\t50\t150\t+\t'
                                      'D\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    culled_alignments = verticall.alignment.cull_redundant_alignments([a, b], 0)
    assert len(culled_alignments) == 2


def test_cull_redundant_alignments_2():
    # Far apart, so no culling.
    a = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    b = verticall.alignment.Alignment('A\t1000\t450\t550\t+\t'
                                      'C\t1000\t450\t550\t100\t100\tAS:i:100\tcg:Z:100=')
    culled_alignments = verticall.alignment.cull_redundant_alignments([a, b], 0)
    assert len(culled_alignments) == 2


def test_cull_redundant_alignments_3():
    # Alignment b is culled because it's contained in alignment a in the query.
    a = verticall.alignment.Alignment('A\t1000\t0\t200\t+\t'
                                      'B\t1000\t0\t200\t100\t100\tAS:i:200\tcg:Z:200=')
    b = verticall.alignment.Alignment('A\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    culled_alignments = verticall.alignment.cull_redundant_alignments([a, b], 0)
    assert len(culled_alignments) == 1
    assert culled_alignments[0].target_name == 'B'


def test_cull_redundant_alignments_4():
    # Alignment b is culled because it's contained in alignment a in the target.
    a = verticall.alignment.Alignment('A\t1000\t0\t200\t+\t'
                                      'C\t1000\t0\t200\t100\t100\tAS:i:200\tcg:Z:200=')
    b = verticall.alignment.Alignment('B\t1000\t50\t150\t+\t'
                                      'C\t1000\t50\t150\t100\t100\tAS:i:100\tcg:Z:100=')
    culled_alignments = verticall.alignment.cull_redundant_alignments([a, b], 0)
    assert len(culled_alignments) == 1
    assert culled_alignments[0].query_name == 'A'


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
    v = 1
    h = 2
    a = 3
    assert verticall.alignment.find_ambiguous_runs([h, h, h, h, h, h]) == []
    assert verticall.alignment.find_ambiguous_runs([v, v, v, v, v, v]) == []
    assert verticall.alignment.find_ambiguous_runs([h, h, a, h, h, h]) == [(2, 3)]
    assert verticall.alignment.find_ambiguous_runs([v, v, v, a, a, v]) == [(3, 5)]
    assert verticall.alignment.find_ambiguous_runs([a, a, a, h, a, h]) == [(0, 3), (4, 5)]
    assert verticall.alignment.find_ambiguous_runs([a, v, a, a, a, a]) == [(0, 1), (2, 6)]
    assert verticall.alignment.find_ambiguous_runs([a, a, a, a, a, a]) == [(0, 6)]


def test_remove_ambiguous():
    v = 1
    h = 2
    a = 3
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


def test_set_up_sliding_windows_1():
    # CIGAR:            11=1X10=1X12=1X1=1X12=
    # simplified CIGAR: ===========X==========X============X=X============
    # windows:          ----------          ----------          ----------
    #                        ----------          ----------
    #                             ----------          ----------
    #                                  ----------          ----------
    a = verticall.alignment.Alignment('A\t1000\t50\t100\t+\t'
                                      'C\t1000\t50\t100\t50\t50\tAS:i:50\t'
                                      'cg:Z:11=1X10=1X12=1X1=1X12=')
    assert a.expanded_cigar == '===========X==========X============X=X============'
    assert a.simplified_cigar == '===========X==========X============X=X============'

    a.set_up_sliding_windows(10, 5)
    assert len(a.windows) == len(a.windows_no_overlap) == len(a.window_differences) == 9
    assert a.windows == [(0, 10), (5, 15), (10, 20), (15, 25), (20, 30),
                         (25, 35), (30, 40), (35, 45), (40, 50)]
    assert a.windows_no_overlap == [(0, 7), (7, 12), (12, 17), (17, 22), (22, 27),
                                    (27, 32), (32, 37), (37, 42), (42, 50)]
    assert a.window_differences == [0, 1, 1, 1, 1, 0, 2, 2, 0]
    assert a.get_max_differences() == 2


def test_set_up_sliding_windows_2():
    # Tests a sliding window too big for this alignment.
    a = verticall.alignment.Alignment('A\t1000\t50\t100\t+\t'
                                      'C\t1000\t50\t100\t50\t50\tAS:i:50\t'
                                      'cg:Z:11=1X10=1X12=1X1=1X12=')
    assert a.expanded_cigar == '===========X==========X============X=X============'
    assert a.simplified_cigar == '===========X==========X============X=X============'

    a.set_up_sliding_windows(100, 10)
    assert len(a.windows) == len(a.windows_no_overlap) == len(a.window_differences) == 0
    assert a.get_max_differences() == 0


def test_paint_sliding_windows_1():
    # CIGAR:            11=1X10=1X12=1X1=1X12=
    # simplified CIGAR: ===========X==========X============X=X============
    # windows:          ----------          ----------          ----------
    #                        ----------          ----------
    #                             ----------          ----------
    #                                  ----------          ----------
    a = verticall.alignment.Alignment('A\t1000\t50\t100\t+\t'
                                      'C\t1000\t50\t100\t50\t50\tAS:i:50\t'
                                      'cg:Z:11=1X10=1X12=1X1=1X12=')
    a.set_up_sliding_windows(10, 5)
    assert a.window_differences == [0, 1, 1, 1, 1, 0, 2, 2, 0]
    thresholds = {'very_low': None, 'low': None, 'high': 1, 'very_high': 2}
    a.paint_sliding_windows(thresholds)
    assert a.window_class_with_amb == [1, 1, 1, 1, 1, 1, 3, 3, 1]
    assert a.window_classifications == [1, 1, 1, 1, 1, 1, 1, 1, 1]
    assert a.get_vertical_blocks(include_ambiguous=False) == [(0, 50)]
    assert a.get_vertical_blocks(include_ambiguous=True) == [(0, 32), (42, 50)]
    assert a.get_horizontal_blocks(include_ambiguous=False) == []
    assert a.get_horizontal_blocks(include_ambiguous=True) == []
    assert a.get_ambiguous_blocks(include_ambiguous=False) == []
    assert a.get_ambiguous_blocks(include_ambiguous=True) == [(32, 42)]
    assert a.get_all_vertical_distances() == [0, 1, 1, 1, 1, 0, 2, 2, 0]
    assert a.get_all_horizontal_distances() == []


def test_paint_sliding_windows_2():
    # CIGAR:            11=1X10=1X12=3X12=
    # simplified CIGAR: ===========X==========X============XXX============
    # windows:          ----------          ----------          ----------
    #                        ----------          ----------
    #                             ----------          ----------
    #                                  ----------          ----------
    a = verticall.alignment.Alignment('A\t1000\t50\t100\t+\t'
                                      'C\t1000\t50\t100\t50\t50\tAS:i:50\t'
                                      'cg:Z:11=1X10=1X12=3X12=')
    a.set_up_sliding_windows(10, 5)
    assert a.window_differences == [0, 1, 1, 1, 1, 0, 3, 3, 0]
    thresholds = {'very_low': None, 'low': None, 'high': 1, 'very_high': 2}
    a.paint_sliding_windows(thresholds)
    assert a.window_class_with_amb == [1, 1, 1, 1, 1, 1, 2, 2, 1]
    assert a.window_classifications == [1, 1, 1, 1, 1, 1, 2, 2, 1]
    assert a.get_vertical_blocks(include_ambiguous=False) == [(0, 32), (42, 50)]
    assert a.get_vertical_blocks(include_ambiguous=True) == [(0, 32), (42, 50)]
    assert a.get_horizontal_blocks(include_ambiguous=False) == [(32, 42)]
    assert a.get_horizontal_blocks(include_ambiguous=True) == [(32, 42)]
    assert a.get_ambiguous_blocks(include_ambiguous=False) == []
    assert a.get_ambiguous_blocks(include_ambiguous=True) == []
    assert a.get_all_vertical_distances() == [0, 1, 1, 1, 1, 0, 0]
    assert a.get_all_horizontal_distances() == [3, 3]


def test_paint_sliding_windows_3():
    # CIGAR:            2=1X2=1X2=1X2=1X2=1X2=1X2=1X2=1X2=1X2=1X11=1X2=1X2=1X2=
    # simplified CIGAR: ==X==X==X==X==X==X==X==X==X==X===========X==X==X==
    # windows:          ----------          ----------          ----------
    #                        ----------          ----------
    #                             ----------          ----------
    #                                  ----------          ----------
    a = verticall.alignment.Alignment('A\t1000\t50\t100\t+\t'
                                      'C\t1000\t50\t100\t50\t50\tAS:i:50\t'
                                      'cg:Z:'
                                      '2=1X2=1X2=1X2=1X2=1X2=1X2=1X2=1X2=1X2=1X11=1X2=1X2=1X2=')
    a.set_up_sliding_windows(10, 5)
    assert a.window_differences == [3, 4, 3, 3, 4, 2, 0, 2, 3]
    thresholds = {'very_low': 1, 'low': 3, 'high': None, 'very_high': None}
    a.paint_sliding_windows(thresholds)
    assert a.window_class_with_amb == [1, 1, 1, 1, 1, 3, 2, 3, 1]
    assert a.window_classifications == [1, 1, 1, 1, 1, 2, 2, 2, 1]
    assert a.get_vertical_blocks(include_ambiguous=False) == [(0, 27), (42, 50)]
    assert a.get_vertical_blocks(include_ambiguous=True) == [(0, 27), (42, 50)]
    assert a.get_horizontal_blocks(include_ambiguous=False) == [(27, 42)]
    assert a.get_horizontal_blocks(include_ambiguous=True) == [(32, 37)]
    assert a.get_ambiguous_blocks(include_ambiguous=False) == []
    assert a.get_ambiguous_blocks(include_ambiguous=True) == [(27, 32), (37, 42)]
    assert a.get_all_vertical_distances() == [3, 4, 3, 3, 4, 3]
    assert a.get_all_horizontal_distances() == [2, 0, 2]
