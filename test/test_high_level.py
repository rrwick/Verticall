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
import re
import tempfile

import verticall.alignment
import verticall.pairwise


def get_random_base():
    return {0: 'A', 1: 'C', 2: 'G', 3: 'T'}[random.randint(0, 3)]


def get_random_different_base(b):
    random_base = get_random_base()
    while b == random_base:
        random_base = get_random_base()
    return random_base


def get_random_seq(seq_len):
    return ''.join(get_random_base() for _ in range(seq_len))


def mutate_seq(seq, divergence):
    change_count = int(round(len(seq) * divergence))
    change_positions = random.sample(range(len(seq)), change_count)
    seq = [b for b in seq]
    for i in change_positions:
        seq[i] = get_random_different_base(seq[i])
    return ''.join(seq)


def reverse_complement(seq):
    rev_bases = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join([rev_bases[x] for x in seq][::-1])


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


def set_up(temp_dir, args, seq_a, seq_b, threads=1):
    temp_dir = pathlib.Path(temp_dir)
    assembly_a = temp_dir / 'A.fasta'
    assembly_b = temp_dir / 'B.fasta'
    write_fasta(assembly_a, 'a', seq_a)
    write_fasta(assembly_b, 'b', seq_b)
    assemblies = verticall.pairwise.find_assemblies(temp_dir)
    verticall.alignment.build_indices(args, assemblies, threads)
    return assembly_a, assembly_b


def assert_regions_approximately_equal(r1, r2, tolerance=0.01):
    """
    Tests if two Verticall region strings contain approximately the same numbers.
    """
    pattern = re.compile(r'[0-9]+')
    numbers_1 = [int(x) for x in re.findall(pattern, r1)]
    numbers_2 = [int(x) for x in re.findall(pattern, r2)]
    assert len(numbers_1) == len(numbers_2)
    for n_1, n_2 in zip(numbers_1, numbers_2):
        assert n_1 == pytest.approx(n_2, tolerance)


def test_no_alignments():
    # Two completely different sequences 500 kbp in length.
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
    # Two identical sequences 500 kbp in length.
    random.seed(0)
    seq_a = get_random_seq(500000)
    seq_b = seq_a
    with tempfile.TemporaryDirectory() as temp_dir:
        args = get_args(temp_dir)
        assembly_a, assembly_b = set_up(temp_dir, args, seq_a, seq_b, 2)
        log_text, table_lines = \
            verticall.pairwise.process_one_pair((args, 'A', 'B', assembly_a, assembly_b))
    assert len(table_lines) == 1
    results = dict(zip(verticall.pairwise.get_table_header().rstrip('\n').split('\t'),
                       table_lines[0].rstrip('\n').split('\t')))

    assert results['assembly_a'] == 'A'
    assert results['assembly_b'] == 'B'
    assert int(results['alignment_count']) == 1
    assert int(results['n50_alignment_length']) == 500000
    assert float(results['aligned_fraction']) == pytest.approx(1.0)
    assert results['result_level'] == 'primary'
    assert float(results['peak_mass']) == pytest.approx(1.0)
    assert float(results['mean_vertical_distance']) == pytest.approx(0.0)
    assert results['r/m'] == 'undef'
    assert results['assembly_a_vertical_regions'] == 'a:0-500000'
    assert results['assembly_a_horizontal_regions'] == ''
    assert results['assembly_a_unaligned_regions'] == ''
    assert results['assembly_b_vertical_regions'] == 'b:0-500000'
    assert results['assembly_b_horizontal_regions'] == ''
    assert results['assembly_b_unaligned_regions'] == ''


def test_identical_different_start_opposite_strand():
    # Two identical sequences 500 kbp in length, but with different start positions and on opposite
    # strands.
    random.seed(0)
    seq_a = get_random_seq(500000)
    seq_b = reverse_complement(seq_a[250000:] + seq_a[:250000])
    with tempfile.TemporaryDirectory() as temp_dir:
        args = get_args(temp_dir)
        assembly_a, assembly_b = set_up(temp_dir, args, seq_a, seq_b, 4)
        log_text, table_lines = \
            verticall.pairwise.process_one_pair((args, 'A', 'B', assembly_a, assembly_b))
    assert len(table_lines) == 1
    results = dict(zip(verticall.pairwise.get_table_header().rstrip('\n').split('\t'),
                       table_lines[0].rstrip('\n').split('\t')))

    assert results['assembly_a'] == 'A'
    assert results['assembly_b'] == 'B'
    assert int(results['alignment_count']) == 2
    assert int(results['n50_alignment_length']) == 250000
    assert float(results['aligned_fraction']) == pytest.approx(1.0)
    assert results['result_level'] == 'primary'
    assert float(results['peak_mass']) == pytest.approx(1.0)
    assert float(results['mean_vertical_distance']) == pytest.approx(0.0)
    assert results['r/m'] == 'undef'
    assert results['assembly_a_vertical_regions'] == 'a:0-500000'
    assert results['assembly_a_horizontal_regions'] == ''
    assert results['assembly_a_unaligned_regions'] == ''
    assert results['assembly_b_vertical_regions'] == 'b:0-500000'
    assert results['assembly_b_horizontal_regions'] == ''
    assert results['assembly_b_unaligned_regions'] == ''


def test_one_unaligned_block():
    # Two identical sequences 500 kbp in length, but with 50 kbp that is completely different.
    # The two sequences are on the same strand.
    random.seed(0)
    seq_a = get_random_seq(500000)
    seq_b = seq_a[0:100000] + get_random_seq(50000) + seq_a[150000:]
    with tempfile.TemporaryDirectory() as temp_dir:
        args = get_args(temp_dir)
        assembly_a, assembly_b = set_up(temp_dir, args, seq_a, seq_b)
        log_text, table_lines = \
            verticall.pairwise.process_one_pair((args, 'A', 'B', assembly_a, assembly_b))
    assert len(table_lines) == 1
    results = dict(zip(verticall.pairwise.get_table_header().rstrip('\n').split('\t'),
                       table_lines[0].rstrip('\n').split('\t')))

    assert results['assembly_a'] == 'A'
    assert results['assembly_b'] == 'B'
    assert int(results['alignment_count']) == 2
    assert int(results['n50_alignment_length']) == pytest.approx(350000, 0.01)
    assert float(results['aligned_fraction']) == pytest.approx(0.9, 0.01)
    assert results['result_level'] == 'primary'
    assert float(results['peak_mass']) == pytest.approx(1.0)
    assert float(results['mean_vertical_distance']) == pytest.approx(0.0)
    assert results['r/m'] == 'undef'
    assert_regions_approximately_equal(results['assembly_a_vertical_regions'],
                                       'a:0-100000,a:150000-500000')
    assert results['assembly_a_horizontal_regions'] == ''
    assert_regions_approximately_equal(results['assembly_a_unaligned_regions'], 'a:100000-150000')
    assert_regions_approximately_equal(results['assembly_b_vertical_regions'],
                                       'b:0-100000,b:150000-500000')
    assert results['assembly_b_horizontal_regions'] == ''
    assert_regions_approximately_equal(results['assembly_b_unaligned_regions'], 'b:100000-150000')


def test_one_unaligned_block_opposite_strand():
    # Two identical sequences 500 kbp in length, but with 50 kbp that is completely different.
    # The two sequences are on opposite strands.
    random.seed(0)
    seq_a = get_random_seq(500000)
    seq_b = reverse_complement(seq_a[0:100000] + get_random_seq(50000) + seq_a[150000:])
    with tempfile.TemporaryDirectory() as temp_dir:
        args = get_args(temp_dir)
        assembly_a, assembly_b = set_up(temp_dir, args, seq_a, seq_b)
        log_text, table_lines = \
            verticall.pairwise.process_one_pair((args, 'A', 'B', assembly_a, assembly_b))
    assert len(table_lines) == 1
    results = dict(zip(verticall.pairwise.get_table_header().rstrip('\n').split('\t'),
                       table_lines[0].rstrip('\n').split('\t')))

    assert results['assembly_a'] == 'A'
    assert results['assembly_b'] == 'B'
    assert int(results['alignment_count']) == 2
    assert int(results['n50_alignment_length']) == pytest.approx(350000, 0.01)
    assert float(results['aligned_fraction']) == pytest.approx(0.9, 0.01)
    assert results['result_level'] == 'primary'
    assert float(results['peak_mass']) == pytest.approx(1.0)
    assert float(results['mean_vertical_distance']) == pytest.approx(0.0)
    assert results['r/m'] == 'undef'
    assert_regions_approximately_equal(results['assembly_a_vertical_regions'],
                                       'a:0-100000,a:150000-500000')
    assert results['assembly_a_horizontal_regions'] == ''
    assert_regions_approximately_equal(results['assembly_a_unaligned_regions'], 'a:100000-150000')
    assert_regions_approximately_equal(results['assembly_b_vertical_regions'],
                                       'b:0-350000,b:400000-500000')
    assert results['assembly_b_horizontal_regions'] == ''
    assert_regions_approximately_equal(results['assembly_b_unaligned_regions'], 'b:350000-400000')


def test_one_high_horizontal_block():
    # Two identical sequences 500 kbp in length, but with 100 kbp that is more divergent.
    # The two sequences are on the same strand.
    random.seed(0)
    seq_a = get_random_seq(500000)
    seq_b = seq_a[0:100000] + mutate_seq(seq_a[100000:200000], 0.01) + seq_a[200000:]
    with tempfile.TemporaryDirectory() as temp_dir:
        args = get_args(temp_dir)
        assembly_a, assembly_b = set_up(temp_dir, args, seq_a, seq_b)
        log_text, table_lines = \
            verticall.pairwise.process_one_pair((args, 'A', 'B', assembly_a, assembly_b))
    assert len(table_lines) == 1
    results = dict(zip(verticall.pairwise.get_table_header().rstrip('\n').split('\t'),
                       table_lines[0].rstrip('\n').split('\t')))

    assert results['assembly_a'] == 'A'
    assert results['assembly_b'] == 'B'
    assert int(results['alignment_count']) == 1
    assert int(results['n50_alignment_length']) == 500000
    assert float(results['aligned_fraction']) == pytest.approx(1.0)
    assert results['result_level'] == 'primary'
    assert float(results['peak_mass']) == pytest.approx(0.8, 0.01)
    assert float(results['mean_vertical_distance']) == pytest.approx(0.0)
    assert results['r/m'] == 'inf'
    assert_regions_approximately_equal(results['assembly_a_vertical_regions'],
                                       'a:0-100000,a:200000-500000')
    assert_regions_approximately_equal(results['assembly_a_horizontal_regions'], 'a:100000-200000')
    assert results['assembly_a_unaligned_regions'] == ''
    assert_regions_approximately_equal(results['assembly_b_vertical_regions'],
                                       'b:0-100000,b:200000-500000')
    assert_regions_approximately_equal(results['assembly_b_horizontal_regions'], 'b:100000-200000')
    assert results['assembly_b_unaligned_regions'] == ''


def test_one_high_horizontal_block_opposite_strand():
    # Two identical sequences 500 kbp in length, but with 100 kbp that is more divergent.
    # The two sequences are on opposite strands.
    random.seed(0)
    seq_a = get_random_seq(500000)
    seq_b = reverse_complement(seq_a[0:100000] + mutate_seq(seq_a[100000:200000], 0.01) +
                               seq_a[200000:])
    with tempfile.TemporaryDirectory() as temp_dir:
        args = get_args(temp_dir)
        assembly_a, assembly_b = set_up(temp_dir, args, seq_a, seq_b)
        log_text, table_lines = \
            verticall.pairwise.process_one_pair((args, 'A', 'B', assembly_a, assembly_b))
    assert len(table_lines) == 1
    results = dict(zip(verticall.pairwise.get_table_header().rstrip('\n').split('\t'),
                       table_lines[0].rstrip('\n').split('\t')))

    assert results['assembly_a'] == 'A'
    assert results['assembly_b'] == 'B'
    assert int(results['alignment_count']) == 1
    assert int(results['n50_alignment_length']) == 500000
    assert float(results['aligned_fraction']) == pytest.approx(1.0)
    assert results['result_level'] == 'primary'
    assert float(results['peak_mass']) == pytest.approx(0.8, 0.01)
    assert float(results['mean_vertical_distance']) == pytest.approx(0.0)
    assert results['r/m'] == 'inf'
    assert_regions_approximately_equal(results['assembly_a_vertical_regions'],
                                       'a:0-100000,a:200000-500000')
    assert_regions_approximately_equal(results['assembly_a_horizontal_regions'], 'a:100000-200000')
    assert results['assembly_a_unaligned_regions'] == ''
    assert_regions_approximately_equal(results['assembly_b_vertical_regions'],
                                       'b:0-300000,b:400000-500000')
    assert_regions_approximately_equal(results['assembly_b_horizontal_regions'], 'b:300000-400000')
    assert results['assembly_b_unaligned_regions'] == ''


def test_one_low_horizontal_block():
    # Two sequences 500 kbp in length that are ~1% divergent, but with 25 kbp that is identical.
    # The two sequences are on the same strand.
    random.seed(0)
    seq_a = get_random_seq(500000)
    seq_b = mutate_seq(seq_a[0:200000], 0.01) + seq_a[200000:225000] + \
        mutate_seq(seq_a[225000:], 0.01)
    with tempfile.TemporaryDirectory() as temp_dir:
        args = get_args(temp_dir)
        assembly_a, assembly_b = set_up(temp_dir, args, seq_a, seq_b)
        log_text, table_lines = \
            verticall.pairwise.process_one_pair((args, 'A', 'B', assembly_a, assembly_b))
    assert len(table_lines) == 1
    results = dict(zip(verticall.pairwise.get_table_header().rstrip('\n').split('\t'),
                       table_lines[0].rstrip('\n').split('\t')))

    assert results['assembly_a'] == 'A'
    assert results['assembly_b'] == 'B'
    assert int(results['alignment_count']) == 1
    assert int(results['n50_alignment_length']) == 500000
    assert float(results['aligned_fraction']) == pytest.approx(1.0)
    assert results['result_level'] == 'primary'
    assert float(results['peak_mass']) == pytest.approx(0.95, 0.01)
    assert float(results['mean_vertical_distance']) == pytest.approx(0.01, 0.01)
    assert float(results['r/m']) == 0.0
    assert_regions_approximately_equal(results['assembly_a_vertical_regions'],
                                       'a:0-200000,a:225000-500000')
    assert_regions_approximately_equal(results['assembly_a_horizontal_regions'], 'a:200000-225000')
    assert results['assembly_a_unaligned_regions'] == ''
    assert_regions_approximately_equal(results['assembly_b_vertical_regions'],
                                       'b:0-200000,b:225000-500000')
    assert_regions_approximately_equal(results['assembly_b_horizontal_regions'], 'b:200000-225000')
    assert results['assembly_b_unaligned_regions'] == ''


def test_one_low_horizontal_block_opposite_strand():
    # Two sequences 500 kbp in length that are ~1% divergent, but with 25 kbp that is identical.
    # The two sequences are on opposite strands.
    random.seed(0)
    seq_a = get_random_seq(500000)
    seq_b = reverse_complement(mutate_seq(seq_a[0:200000], 0.01) + seq_a[200000:225000] +
                               mutate_seq(seq_a[225000:], 0.01))
    with tempfile.TemporaryDirectory() as temp_dir:
        args = get_args(temp_dir)
        assembly_a, assembly_b = set_up(temp_dir, args, seq_a, seq_b)
        log_text, table_lines = \
            verticall.pairwise.process_one_pair((args, 'A', 'B', assembly_a, assembly_b))
    assert len(table_lines) == 1
    results = dict(zip(verticall.pairwise.get_table_header().rstrip('\n').split('\t'),
                       table_lines[0].rstrip('\n').split('\t')))

    assert results['assembly_a'] == 'A'
    assert results['assembly_b'] == 'B'
    assert int(results['alignment_count']) == 1
    assert int(results['n50_alignment_length']) == 500000
    assert float(results['aligned_fraction']) == pytest.approx(1.0)
    assert results['result_level'] == 'primary'
    assert float(results['peak_mass']) == pytest.approx(0.95, 0.01)
    assert float(results['mean_vertical_distance']) == pytest.approx(0.01, 0.01)
    assert float(results['r/m']) == 0.0
    assert_regions_approximately_equal(results['assembly_a_vertical_regions'],
                                       'a:0-200000,a:225000-500000')
    assert_regions_approximately_equal(results['assembly_a_horizontal_regions'], 'a:200000-225000')
    assert results['assembly_a_unaligned_regions'] == ''
    assert_regions_approximately_equal(results['assembly_b_vertical_regions'],
                                       'b:0-275000,b:300000-500000')
    assert_regions_approximately_equal(results['assembly_b_horizontal_regions'], 'b:275000-300000')
    assert results['assembly_b_unaligned_regions'] == ''
