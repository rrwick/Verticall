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

import gzip
import pytest

import verticall.misc


def test_get_compression_type_1():
    assert verticall.misc.get_compression_type('test/test_misc/test.txt') == 'plain'


def test_get_compression_type_2():
    assert verticall.misc.get_compression_type('test/test_misc/test.gz') == 'gz'


def test_get_compression_type_3():
    with pytest.raises(SystemExit) as e:
        verticall.misc.get_compression_type('test/test_misc/test.bz2')
    assert 'cannot use bzip2' in str(e.value)


def test_get_compression_type_4():
    with pytest.raises(SystemExit) as e:
        verticall.misc.get_compression_type('test/test_misc/test.zip')
    assert 'cannot use zip' in str(e.value)


def test_get_open_func_1():
    assert verticall.misc.get_open_func('test/test_misc/test.txt') == open


def test_get_open_func_2():
    assert verticall.misc.get_open_func('test/test_misc/test.gz') == gzip.open


def test_get_sequence_file_type_01():
    assert verticall.misc.get_sequence_file_type('test/test_misc/test.fasta') == 'FASTA'


def test_get_sequence_file_type_02():
    assert verticall.misc.get_sequence_file_type('test/test_misc/test.fastq') == 'FASTQ'


def test_get_sequence_file_type_03():
    assert verticall.misc.get_sequence_file_type('test/test_misc/test.gfa') == 'GFA'


def test_get_sequence_file_type_04():
    assert verticall.misc.get_sequence_file_type('test/test_misc/test.txt') == 'unknown'


def test_get_sequence_file_type_05():
    assert verticall.misc.get_sequence_file_type('test/test_misc/empty') == 'unknown'


def test_get_sequence_file_type_06():
    assert verticall.misc.get_sequence_file_type('test/test_misc/test.fasta.gz') == 'FASTA'


def test_get_sequence_file_type_07():
    assert verticall.misc.get_sequence_file_type('test/test_misc/test.fastq.gz') == 'FASTQ'


def test_get_sequence_file_type_08():
    assert verticall.misc.get_sequence_file_type('test/test_misc/test.fastq.gz') == 'FASTQ'


def test_get_sequence_file_type_09():
    assert verticall.misc.get_sequence_file_type('test/test_misc/test.gfa.gz') == 'GFA'


def test_get_sequence_file_type_10():
    assert verticall.misc.get_sequence_file_type('test/test_misc/not_unicode') == 'unknown'


def test_iterate_fasta_1():
    fasta = list(verticall.misc.iterate_fasta('test/test_misc/test.fasta'))
    assert fasta == [('A', 'TTGCCTGTAGTCGGGACCCC'), ('B', 'ATTCTCAGAATGGCGTAGTA'),
                     ('C', 'TACGCAGCTACG')]
    fasta = list(verticall.misc.iterate_fasta('test/test_misc/test.fasta', include_info=True))
    assert fasta == [('A', 'info', 'TTGCCTGTAGTCGGGACCCC'),
                     ('B', 'stuff and more stuff', 'ATTCTCAGAATGGCGTAGTA'),
                     ('C', '', 'TACGCAGCTACG')]


def test_iterate_fasta_2():
    fasta = list(verticall.misc.iterate_fasta('test/test_misc/test.fasta.gz'))
    assert fasta == [('A', 'TTGCCTGTAGTCGGGACCCC'), ('B', 'ATTCTCAGAATGGCGTAGTA'),
                     ('C', 'TACGCAGCTACG')]
    fasta = list(verticall.misc.iterate_fasta('test/test_misc/test.fasta.gz', include_info=True))
    assert fasta == [('A', 'info', 'TTGCCTGTAGTCGGGACCCC'),
                     ('B', 'stuff and more stuff', 'ATTCTCAGAATGGCGTAGTA'),
                     ('C', '', 'TACGCAGCTACG')]


def test_get_default_thread_count():
    assert 1 <= verticall.misc.get_default_thread_count() <= 16


def test_get_n50():
    assert verticall.misc.get_n50([1, 2, 3, 4, 1000]) == 1000
    assert verticall.misc.get_n50([12, 23455, 15, 12433, 15343, 9, 10]) == 15343
    assert verticall.misc.get_n50([]) == 0


def test_get_window_count():
    # ----------      ----------      ----------
    #         ----------      ----------      ----------
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    assert verticall.misc.get_window_count(50, 10, 8) == 6

    # ----------       ----------       ----------
    #          ----------       ----------
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    assert verticall.misc.get_window_count(50, 10, 9) == 5

    # ----------          ----------          ----------
    #           ----------          ----------
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    assert verticall.misc.get_window_count(50, 10, 10) == 5

    # ----------            ----------
    #            ----------            ----------
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    assert verticall.misc.get_window_count(50, 10, 11) == 4


def test_get_window_coverage():
    assert verticall.misc.get_window_coverage(10, 8, 5) == 42


def test_split_list():
    lst = [0, 1, 2, 3, 4, 5, 6, 7]
    assert verticall.misc.split_list(lst, 1) == [[0, 1, 2, 3, 4, 5, 6, 7]]
    assert verticall.misc.split_list(lst, 2) == [[0, 1, 2, 3], [4, 5, 6, 7]]
    assert verticall.misc.split_list(lst, 3) == [[0, 1, 2], [3, 4, 5], [6, 7]]
    assert verticall.misc.split_list(lst, 4) == [[0, 1], [2, 3], [4, 5], [6, 7]]
    assert verticall.misc.split_list(lst, 8) == [[0], [1], [2], [3], [4], [5], [6], [7]]


def test_contains_ambiguous_bases():
    assert not verticall.misc.contains_ambiguous_bases('ACGATCGACTACG')
    assert not verticall.misc.contains_ambiguous_bases('acgatcgacgac')
    assert not verticall.misc.contains_ambiguous_bases('AAAAcccccTTTTggg')
    assert verticall.misc.contains_ambiguous_bases('ACGACTAGCNACTAGCACT')
    assert verticall.misc.contains_ambiguous_bases('ACGACTAGCMACTAGCACT')
    assert verticall.misc.contains_ambiguous_bases('ACGACTAGCRACTAGCACT')
    assert verticall.misc.contains_ambiguous_bases('NNNNNNNN')
    assert verticall.misc.contains_ambiguous_bases('DSHIFUSDJFSDOFJ')


def test_list_differences():
    a = [1, 2, 3]
    b = [1, 2, 3]
    in_both, in_a_not_b, in_b_not_a = verticall.misc.list_differences(a, b)
    assert in_both == [1, 2, 3]
    assert in_a_not_b == []
    assert in_b_not_a == []

    a = [1, 2, 3, 4, 5]
    b = [1, 2, 3]
    in_both, in_a_not_b, in_b_not_a = verticall.misc.list_differences(a, b)
    assert in_both == [1, 2, 3]
    assert in_a_not_b == [4, 5]
    assert in_b_not_a == []

    a = [1, 2, 3]
    b = [1, 2, 3, 4, 5]
    in_both, in_a_not_b, in_b_not_a = verticall.misc.list_differences(a, b)
    assert in_both == [1, 2, 3]
    assert in_a_not_b == []
    assert in_b_not_a == [4, 5]

    a = [1, 2, 3]
    b = [4, 5, 6]
    in_both, in_a_not_b, in_b_not_a = verticall.misc.list_differences(a, b)
    assert in_both == []
    assert in_a_not_b == [1, 2, 3]
    assert in_b_not_a == [4, 5, 6]
