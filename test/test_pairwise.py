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
import tempfile

import verticall.pairwise


def test_welcome_message(capsys):
    Args = collections.namedtuple('Args', ['reference', 'part'])

    verticall.pairwise.welcome_message(Args(reference=None, part='1/1'))
    _, err = capsys.readouterr()
    assert 'Verticall pairwise' in err
    assert '--part' not in err

    verticall.pairwise.welcome_message(Args(reference=None, part='1/2'))
    _, err = capsys.readouterr()
    assert 'Verticall pairwise' in err
    assert '--part' in err


def test_finished_message(capsys):
    verticall.pairwise.finished_message(True)
    _, err = capsys.readouterr()
    assert 'Finished' in err


def test_parse_part():
    assert verticall.pairwise.parse_part('1/1') == (0, 1)
    assert verticall.pairwise.parse_part('1/2') == (0, 2)
    assert verticall.pairwise.parse_part('2/2') == (1, 2)
    assert verticall.pairwise.parse_part('1/10') == (0, 10)
    assert verticall.pairwise.parse_part('5/10') == (4, 10)
    assert verticall.pairwise.parse_part('10/10') == (9, 10)

    with pytest.raises(SystemExit) as e:
        verticall.pairwise.parse_part('0/1')
    assert 'the numerator of --part must be a positive integer' in str(e.value)

    with pytest.raises(SystemExit) as e:
        verticall.pairwise.parse_part('1/0')
    assert 'the denominator of --part must be a positive integer' in str(e.value)

    with pytest.raises(SystemExit) as e:
        verticall.pairwise.parse_part('2/1')
    assert 'the numerator of --part must be less than or equal to the denominator' in str(e.value)


def test_get_arg_list_1():
    """
    Make sure that we get the same arguments regardless of how we split them into parts.
    """
    get_arg_list = verticall.pairwise.get_arg_list
    Args = collections.namedtuple('Args', ['part', 'existing_tsv'])
    assemblies = [('a', 'a.fasta'), ('b', 'b.fasta'), ('c', 'c.fasta'), ('d', 'd.fasta'),
                  ('e', 'e.fasta'), ('f', 'f.fasta'), ('g', 'g.fasta'), ('h', 'h.fasta')]

    arg_list_1 = get_arg_list(Args(part='1/1', existing_tsv=None), assemblies, None)
    assert len(arg_list_1) == 56

    arg_list_1_2 = get_arg_list(Args(part='1/2', existing_tsv=None), assemblies, None)
    arg_list_2_2 = get_arg_list(Args(part='2/2', existing_tsv=None), assemblies, None)
    arg_list_2 = arg_list_1_2 + arg_list_2_2
    assert [a[1:] for a in arg_list_1] == [a[1:] for a in arg_list_2]

    arg_list_1_3 = get_arg_list(Args(part='1/3', existing_tsv=None), assemblies, None)
    arg_list_2_3 = get_arg_list(Args(part='2/3', existing_tsv=None), assemblies, None)
    arg_list_3_3 = get_arg_list(Args(part='3/3', existing_tsv=None), assemblies, None)
    arg_list_3 = arg_list_1_3 + arg_list_2_3 + arg_list_3_3
    assert [a[1:] for a in arg_list_1] == [a[1:] for a in arg_list_3]

    arg_list_1_4 = get_arg_list(Args(part='1/4', existing_tsv=None), assemblies, None)
    arg_list_2_4 = get_arg_list(Args(part='2/4', existing_tsv=None), assemblies, None)
    arg_list_3_4 = get_arg_list(Args(part='3/4', existing_tsv=None), assemblies, None)
    arg_list_4_4 = get_arg_list(Args(part='4/4', existing_tsv=None), assemblies, None)
    arg_list_4 = arg_list_1_4 + arg_list_2_4 + arg_list_3_4 + arg_list_4_4
    assert [a[1:] for a in arg_list_1] == [a[1:] for a in arg_list_4]

    arg_list_1_5 = get_arg_list(Args(part='1/5', existing_tsv=None), assemblies, None)
    arg_list_2_5 = get_arg_list(Args(part='2/5', existing_tsv=None), assemblies, None)
    arg_list_3_5 = get_arg_list(Args(part='3/5', existing_tsv=None), assemblies, None)
    arg_list_4_5 = get_arg_list(Args(part='4/5', existing_tsv=None), assemblies, None)
    arg_list_5_5 = get_arg_list(Args(part='5/5', existing_tsv=None), assemblies, None)
    arg_list_5 = arg_list_1_5 + arg_list_2_5 + arg_list_3_5 + arg_list_4_5 + arg_list_5_5
    assert [a[1:] for a in arg_list_1] == [a[1:] for a in arg_list_5]


def test_get_arg_list_2():
    Args = collections.namedtuple('Args', ['part', 'existing_tsv'])
    assemblies = [('a', 'a.fasta'), ('b', 'b.fasta'), ('c', 'c.fasta'), ('d', 'd.fasta'),
                  ('e', 'e.fasta'), ('f', 'f.fasta'), ('g', 'g.fasta'), ('h', 'h.fasta')]
    reference = ('ref', 'ref.fasta')

    arg_list = verticall.pairwise.get_arg_list(Args(part='1/1', existing_tsv=None), assemblies,
                                               reference)
    assert len(arg_list) == 8


def test_get_arg_list_3():
    Args = collections.namedtuple('Args', ['part', 'existing_tsv'])
    assemblies = [('a', 'a.fasta'), ('b', 'b.fasta'), ('c', 'c.fasta'), ('d', 'd.fasta'),
                  ('e', 'e.fasta'), ('f', 'f.fasta'), ('g', 'g.fasta'), ('h', 'h.fasta')]

    arg_list = verticall.pairwise.get_arg_list(Args(part='1/1', existing_tsv=None), assemblies,
                                               None)
    assert len(arg_list) == 56


def test_get_arg_list_4():
    Args = collections.namedtuple('Args', ['part', 'existing_tsv'])
    assemblies = [('a', 'a.fasta'), ('b', 'b.fasta'), ('c', 'c.fasta'), ('d', 'd.fasta'),
                  ('e', 'e.fasta'), ('f', 'f.fasta'), ('g', 'g.fasta'), ('h', 'h.fasta')]

    with tempfile.TemporaryDirectory() as temp_dir:
        out_file = pathlib.Path(temp_dir) / 'existing.tsv'
        verticall.pairwise.create_output_dir_if_needed(out_file)
        with open(out_file, 'wt') as f:
            f.write('assembly_a\tassembly_b\n')
            f.write('a\tb\n')
            f.write('a\tc\n')
            f.write('a\td\n')

        arg_list = verticall.pairwise.get_arg_list(Args(part='1/1', existing_tsv=out_file),
                                                   assemblies, None)
        assert len(arg_list) == 53


def test_find_assemblies_1():
    assembly_dir = pathlib.Path('test/test_pairwise/assemblies')
    assemblies = verticall.pairwise.find_assemblies(assembly_dir)
    assert [a[0] for a in assemblies] == ['a', 'b', 'c', 'd', 'e', 'f']
    assemblies = verticall.pairwise.find_assemblies(assembly_dir,
                                                    extensions=['.fasta', '.fasta.gz'])
    assert [a[0] for a in assemblies] == ['a', 'b']
    assemblies = verticall.pairwise.find_assemblies(assembly_dir, extensions=['.other'])
    assert [a[0] for a in assemblies] == ['g']


def test_find_assemblies_2():
    assembly_dir = pathlib.Path('test/test_pairwise/assemblies_dup_name')
    with pytest.raises(SystemExit) as e:
        verticall.pairwise.find_assemblies(assembly_dir)
    assert 'duplicate' in str(e.value)


def test_check_assemblies_1():
    assemblies = [('a', pathlib.Path('test/test_pairwise/assemblies/a.fasta'))]
    verticall.pairwise.check_assemblies(assemblies, 1, None)

    assemblies = [('b', pathlib.Path('test/test_pairwise/assemblies/b.fasta.gz'))]
    with pytest.raises(SystemExit) as e:
        verticall.pairwise.check_assemblies(assemblies, 1, None)
    assert 'duplicate contig names' in str(e.value)

    assemblies = [('c', pathlib.Path('test/test_pairwise/assemblies/c.fna'))]
    with pytest.raises(SystemExit) as e:
        verticall.pairwise.check_assemblies(assemblies, 1, None)
    assert 'ambiguous' in str(e.value)


def test_check_assemblies_2():
    assemblies = [('a', pathlib.Path('test/test_pairwise/assemblies/a.fasta'))]
    verticall.pairwise.check_assemblies(assemblies, 4, None)

    assemblies = [('b', pathlib.Path('test/test_pairwise/assemblies/b.fasta.gz'))]
    with pytest.raises(SystemExit) as e:
        verticall.pairwise.check_assemblies(assemblies, 4, None)
    assert 'duplicate contig names' in str(e.value)

    assemblies = [('c', pathlib.Path('test/test_pairwise/assemblies/c.fna'))]
    with pytest.raises(SystemExit) as e:
        verticall.pairwise.check_assemblies(assemblies, 4, None)
    assert 'ambiguous' in str(e.value)


def test_check_assemblies_3():
    assemblies = [('a', pathlib.Path('test/test_pairwise/assemblies/a.fasta')),
                  ('b', pathlib.Path('test/test_pairwise/assemblies/b.fasta.gz')),
                  ('c', pathlib.Path('test/test_pairwise/assemblies/c.fna'))]
    with pytest.raises(SystemExit) as e:
        verticall.pairwise.check_assemblies(assemblies, 1, None)
    assert 'duplicate contig names' in str(e.value)
    assert 'ambiguous' in str(e.value)


def test_check_assemblies_4():
    assemblies = [('a', pathlib.Path('test/test_pairwise/assemblies/a.fasta')),
                  ('b', pathlib.Path('test/test_pairwise/assemblies/b.fasta.gz')),
                  ('c', pathlib.Path('test/test_pairwise/assemblies/c.fna'))]
    with pytest.raises(SystemExit) as e:
        verticall.pairwise.check_assemblies(assemblies, 4, None)
    assert 'duplicate contig names' in str(e.value)
    assert 'ambiguous' in str(e.value)


def test_check_assemblies_5():
    reference = ('a', pathlib.Path('test/test_pairwise/assemblies/a.fasta'))

    assemblies = [('b', pathlib.Path('test/test_pairwise/assemblies/b.fasta.gz'))]
    with pytest.raises(SystemExit) as e:
        verticall.pairwise.check_assemblies(assemblies, 1, reference)
    assert 'duplicate contig names' in str(e.value)

    assemblies = [('c', pathlib.Path('test/test_pairwise/assemblies/c.fna'))]
    with pytest.raises(SystemExit) as e:
        verticall.pairwise.check_assemblies(assemblies, 4, reference)
    assert 'ambiguous' in str(e.value)


def test_create_output_dir_if_needed_1():
    with tempfile.TemporaryDirectory() as temp_dir:
        out_file = pathlib.Path(temp_dir) / 'out.tsv'
        verticall.pairwise.create_output_dir_if_needed(out_file)
        with open(out_file, 'wt') as f:
            f.write('test')


def test_create_output_dir_if_needed_2():
    with tempfile.TemporaryDirectory() as temp_dir:
        out_file = pathlib.Path(temp_dir) / 'dir1' / 'out.tsv'
        verticall.pairwise.create_output_dir_if_needed(out_file)
        with open(out_file, 'wt') as f:
            f.write('test')


def test_create_output_dir_if_needed_3():
    with tempfile.TemporaryDirectory() as temp_dir:
        out_file = pathlib.Path(temp_dir) / 'dir1' / 'dir2' / 'out.tsv'
        verticall.pairwise.create_output_dir_if_needed(out_file)
        with open(out_file, 'wt') as f:
            f.write('test')


def test_load_existing_pairs_1():
    tsv = pathlib.Path('test/test_matrix/pairwise.tsv')
    pairs = verticall.pairwise.load_existing_pairs(tsv)
    assert len(pairs) == 56


def test_load_existing_pairs_2():
    tsv = pathlib.Path('test/test_matrix/pairwise_blank_lines.tsv')
    pairs = verticall.pairwise.load_existing_pairs(tsv)
    assert len(pairs) == 56


def test_load_existing_pairs_3():
    tsv = pathlib.Path('test/test_matrix/empty.tsv')
    pairs = verticall.pairwise.load_existing_pairs(tsv)
    assert len(pairs) == 0


def test_load_existing_pairs_4():
    tsv = pathlib.Path('test/test_matrix/empty_blank_lines.tsv')
    pairs = verticall.pairwise.load_existing_pairs(tsv)
    assert len(pairs) == 0


def test_load_existing_pairs_5():
    tsv = pathlib.Path('test/test_matrix/header_only.tsv')
    pairs = verticall.pairwise.load_existing_pairs(tsv)
    assert len(pairs) == 0
