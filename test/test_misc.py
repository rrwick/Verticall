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

import gzip
import pytest

import phylo.misc


def test_get_compression_type_1():
    assert phylo.misc.get_compression_type('test/test_misc/test.txt') == 'plain'


def test_get_compression_type_2():
    assert phylo.misc.get_compression_type('test/test_misc/test.gz') == 'gz'


def test_get_compression_type_3():
    with pytest.raises(SystemExit) as e:
        phylo.misc.get_compression_type('test/test_misc/test.bz2')
    assert 'cannot use bzip2' in str(e.value)


def test_get_compression_type_4():
    with pytest.raises(SystemExit) as e:
        phylo.misc.get_compression_type('test/test_misc/test.zip')
    assert 'cannot use zip' in str(e.value)


def test_get_open_func_1():
    assert phylo.misc.get_open_func('test/test_misc/test.txt') == open


def test_get_open_func_2():
    assert phylo.misc.get_open_func('test/test_misc/test.gz') == gzip.open


def test_get_sequence_file_type_01():
    assert phylo.misc.get_sequence_file_type('test/test_misc/test.fasta') == 'FASTA'


def test_get_sequence_file_type_02():
    assert phylo.misc.get_sequence_file_type('test/test_misc/test.fastq') == 'FASTQ'


def test_get_sequence_file_type_03():
    assert phylo.misc.get_sequence_file_type('test/test_misc/test.gfa') == 'GFA'


def test_get_sequence_file_type_04():
    assert phylo.misc.get_sequence_file_type('test/test_misc/test.txt') == 'unknown'


def test_get_sequence_file_type_05():
    assert phylo.misc.get_sequence_file_type('test/test_misc/empty') == 'unknown'


def test_get_sequence_file_type_06():
    assert phylo.misc.get_sequence_file_type('test/test_misc/test.fasta.gz') == 'FASTA'


def test_get_sequence_file_type_07():
    assert phylo.misc.get_sequence_file_type('test/test_misc/test.fastq.gz') == 'FASTQ'


def test_get_sequence_file_type_08():
    assert phylo.misc.get_sequence_file_type('test/test_misc/test.fastq.gz') == 'FASTQ'


def test_get_sequence_file_type_09():
    assert phylo.misc.get_sequence_file_type('test/test_misc/test.gfa.gz') == 'GFA'


def test_get_sequence_file_type_10():
    assert phylo.misc.get_sequence_file_type('test/test_misc/not_unicode') == 'unknown'


def test_iterate_fasta_1():
    fasta = list(phylo.misc.iterate_fasta('test/test_misc/test.fasta'))
    assert fasta == [('A', 'TTGCCTGTAGTCGGGACCCC'), ('B', 'ATTCTCAGAATGGCGTAGTA')]


def test_iterate_fasta_2():
    fasta = list(phylo.misc.iterate_fasta('test/test_misc/test.fasta.gz'))
    assert fasta == [('A', 'TTGCCTGTAGTCGGGACCCC'), ('B', 'ATTCTCAGAATGGCGTAGTA')]


def test_check_assembly_file_type_1():
    phylo.misc.check_assembly_file_type('test/test_misc/test.fasta')
    phylo.misc.check_assembly_file_type('test/test_misc/test.fasta.gz')


def test_check_assembly_file_type_2():
    with pytest.raises(SystemExit) as e:
        phylo.misc.check_assembly_file_type('test/test_misc/test.fastq')
    assert 'is not in FASTA format' in str(e.value)
    with pytest.raises(SystemExit) as e:
        phylo.misc.check_assembly_file_type('test/test_misc/test.fastq.gz')
    assert 'is not in FASTA format' in str(e.value)


def test_check_assembly_file_type_3():
    with pytest.raises(SystemExit) as e:
        phylo.misc.check_assembly_file_type('test/test_misc/test.gfa')
    assert 'is not in FASTA format' in str(e.value)
    with pytest.raises(SystemExit) as e:
        phylo.misc.check_assembly_file_type('test/test_misc/test.gfa.gz')
    assert 'is not in FASTA format' in str(e.value)


def test_get_default_thread_count():
    assert 1 <= phylo.misc.get_default_thread_count() <= 16


def test_get_n50():
    assert phylo.misc.get_n50([1, 2, 3, 4, 1000]) == 1000
    assert phylo.misc.get_n50([12, 23455, 15, 12433, 15343, 9, 10]) == 15343
    assert phylo.misc.get_n50([]) == 0


def test_get_window_count():
    # ----------      ----------      ----------
    #         ----------      ----------      ----------
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    assert phylo.misc.get_window_count(50, 10, 8) == 6

    # ----------       ----------       ----------
    #          ----------       ----------
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    assert phylo.misc.get_window_count(50, 10, 9) == 5

    # ----------          ----------          ----------
    #           ----------          ----------
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    assert phylo.misc.get_window_count(50, 10, 10) == 5

    # ----------            ----------
    #            ----------            ----------
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    assert phylo.misc.get_window_count(50, 10, 11) == 4


def test_get_window_coverage():
    assert phylo.misc.get_window_coverage(10, 8, 5) == 42
