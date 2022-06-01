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

import collections
import pathlib
import pytest
import random
import tempfile

import verticall.alignment
import verticall.pairwise


def get_random_base():
    return {0: 'A', 1: 'C', 2: 'G', 3: 'T'}[random.randint(0, 3)]


def get_random_seq(seq_len):
    return ''.join(get_random_base() for _ in range(seq_len))


def write_fasta(filename, name, seq):
    with open(filename, 'wt') as f:
        f.write(f'>{name}\n')
        f.write(f'{seq}\n')


def get_args(in_dir, smoothing_factor=0.8, secondary=0.7, ignore_indels=False, allowed_overlap=100,
             window_count=50000, window_size=None, verbose=False, align_options='-x asm20',
             index_options='-k15 -w10'):
    Args = collections.namedtuple('Args', ['smoothing_factor', 'secondary', 'ignore_indels',
                                           'allowed_overlap', 'window_count', 'window_size',
                                           'in_dir', 'verbose', 'align_options', 'index_options'])
    return Args(smoothing_factor=smoothing_factor, secondary=secondary, ignore_indels=ignore_indels,
                allowed_overlap=allowed_overlap, window_count=window_count, window_size=window_size,
                in_dir=pathlib.Path(in_dir), verbose=verbose, align_options=align_options,
                index_options=index_options)


def get_header_index():
    header_line = verticall.pairwise.get_table_header()
    parts = header_line.strip().split('\t')
    return {p: i for i, p in enumerate(parts)}


def set_up(temp_dir, args, seq_a, seq_b):
    temp_dir = pathlib.Path(temp_dir)
    assembly_a = temp_dir / 'A.fasta'
    assembly_b = temp_dir / 'B.fasta'
    write_fasta(assembly_a, 'a', seq_a)
    write_fasta(assembly_b, 'b', seq_b)
    assemblies = verticall.pairwise.find_assemblies(temp_dir)
    verticall.alignment.build_indices(args, assemblies)
    return assembly_a, assembly_b


def test_no_alignments():
    # Tests two completely different sequences 500 kbp in length.
    random.seed(0)
    seq_a = get_random_seq(500000)
    seq_b = get_random_seq(500000)
    with tempfile.TemporaryDirectory() as temp_dir:
        args = get_args(temp_dir)
        assembly_a, assembly_b = set_up(temp_dir, args, seq_a, seq_b)
        log_text, table_lines = \
            verticall.pairwise.process_one_pair((args, 'A', 'B', assembly_a, assembly_b))
    assert len(table_lines) == 0
    assert 'no alignments found' in '\n'.join(log_text)


def test_identical():
    # Tests two identical sequences 500 kbp in length.
    random.seed(0)
    seq_a = get_random_seq(500000)
    seq_b = seq_a
    with tempfile.TemporaryDirectory() as temp_dir:
        args = get_args(temp_dir)
        assembly_a, assembly_b = set_up(temp_dir, args, seq_a, seq_b)
        log_text, table_lines = \
            verticall.pairwise.process_one_pair((args, 'A', 'B', assembly_a, assembly_b))
    assert len(table_lines) == 1
    results = table_lines[0].rstrip('\n').split('\t')
    header = get_header_index()

    assert results[header['assembly_a']] == 'A'
    assert results[header['assembly_b']] == 'B'
    assert int(results[header['alignment_count']]) == 1
    assert int(results[header['n50_alignment_length']]) == 500000
    assert float(results[header['aligned_fraction']]) == pytest.approx(1.0)
    assert results[header['result_level']] == 'primary'
    assert float(results[header['peak_mass']]) == pytest.approx(1.0)
    assert float(results[header['mean_vertical_distance']]) == pytest.approx(0.0)
    assert float(results[header['r/m']]) == pytest.approx(0.0)
    assert results[header['assembly_a_vertical_regions']] == 'a:0-500000'
    assert results[header['assembly_a_horizontal_regions']] == ''
    assert results[header['assembly_a_unaligned_regions']] == ''
    assert results[header['assembly_b_vertical_regions']] == 'b:0-500000'
    assert results[header['assembly_b_horizontal_regions']] == ''
    assert results[header['assembly_b_unaligned_regions']] == ''
