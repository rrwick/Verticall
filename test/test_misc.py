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


def test_get_generator_1():
    assert phylo.misc.get_generator('test/test_misc/test.fasta') == phylo.misc.iterate_fasta


def test_get_generator_2():
    assert phylo.misc.get_generator('test/test_misc/test.gfa') == phylo.misc.iterate_gfa


def test_get_generator_3():
    assert phylo.misc.get_generator('test/test_misc/test.fasta.gz') == phylo.misc.iterate_fasta


def test_get_generator_4():
    assert phylo.misc.get_generator('test/test_misc/test.gfa.gz') == phylo.misc.iterate_gfa


def test_iterate_fasta_1():
    fasta = list(phylo.misc.iterate_fasta('test/test_misc/test.fasta'))
    assert fasta == [('A', 'TTGCCTGTAGTCGGGACCCC'), ('B', 'ATTCTCAGAATGGCGTAGTA')]


def test_iterate_fasta_2():
    fasta = list(phylo.misc.iterate_fasta('test/test_misc/test.fasta.gz'))
    assert fasta == [('A', 'TTGCCTGTAGTCGGGACCCC'), ('B', 'ATTCTCAGAATGGCGTAGTA')]


def test_iterate_gfa_1():
    gfa = list(phylo.misc.iterate_gfa('test/test_misc/test.gfa'))
    assert gfa == [('1', 'ACGTACGT'), ('2', 'TGCA')]


def test_iterate_gfa_2():
    gfa = list(phylo.misc.iterate_gfa('test/test_misc/test.gfa.gz'))
    assert gfa == [('1', 'ACGTACGT'), ('2', 'TGCA')]
