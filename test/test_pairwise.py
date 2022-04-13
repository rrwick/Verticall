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
    verticall.pairwise.finished_message()
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
    Args = collections.namedtuple('Args', ['part'])
    assemblies = [('a', 'a.fasta'), ('b', 'b.fasta'), ('c', 'c.fasta'), ('d', 'd.fasta'),
                  ('e', 'e.fasta'), ('f', 'f.fasta'), ('g', 'g.fasta'), ('h', 'h.fasta')]

    arg_list_1 = verticall.pairwise.get_arg_list(Args(part='1/1'), assemblies, None)
    assert len(arg_list_1) == 56

    arg_list_1_2 = verticall.pairwise.get_arg_list(Args(part='1/2'), assemblies, None)
    arg_list_2_2 = verticall.pairwise.get_arg_list(Args(part='2/2'), assemblies, None)
    arg_list_2 = arg_list_1_2 + arg_list_2_2
    assert [a[1:] for a in arg_list_1] == [a[1:] for a in arg_list_2]

    arg_list_1_3 = verticall.pairwise.get_arg_list(Args(part='1/3'), assemblies, None)
    arg_list_2_3 = verticall.pairwise.get_arg_list(Args(part='2/3'), assemblies, None)
    arg_list_3_3 = verticall.pairwise.get_arg_list(Args(part='3/3'), assemblies, None)
    arg_list_3 = arg_list_1_3 + arg_list_2_3 + arg_list_3_3
    assert [a[1:] for a in arg_list_1] == [a[1:] for a in arg_list_3]

    arg_list_1_4 = verticall.pairwise.get_arg_list(Args(part='1/4'), assemblies, None)
    arg_list_2_4 = verticall.pairwise.get_arg_list(Args(part='2/4'), assemblies, None)
    arg_list_3_4 = verticall.pairwise.get_arg_list(Args(part='3/4'), assemblies, None)
    arg_list_4_4 = verticall.pairwise.get_arg_list(Args(part='4/4'), assemblies, None)
    arg_list_4 = arg_list_1_4 + arg_list_2_4 + arg_list_3_4 + arg_list_4_4
    assert [a[1:] for a in arg_list_1] == [a[1:] for a in arg_list_4]

    arg_list_1_5 = verticall.pairwise.get_arg_list(Args(part='1/5'), assemblies, None)
    arg_list_2_5 = verticall.pairwise.get_arg_list(Args(part='2/5'), assemblies, None)
    arg_list_3_5 = verticall.pairwise.get_arg_list(Args(part='3/5'), assemblies, None)
    arg_list_4_5 = verticall.pairwise.get_arg_list(Args(part='4/5'), assemblies, None)
    arg_list_5_5 = verticall.pairwise.get_arg_list(Args(part='5/5'), assemblies, None)
    arg_list_5 = arg_list_1_5 + arg_list_2_5 + arg_list_3_5 + arg_list_4_5 + arg_list_5_5
    assert [a[1:] for a in arg_list_1] == [a[1:] for a in arg_list_5]


def test_get_arg_list_2():
    Args = collections.namedtuple('Args', ['part'])
    assemblies = [('a', 'a.fasta'), ('b', 'b.fasta'), ('c', 'c.fasta'), ('d', 'd.fasta'),
                  ('e', 'e.fasta'), ('f', 'f.fasta'), ('g', 'g.fasta'), ('h', 'h.fasta')]
    reference = ('ref', 'ref.fasta')

    arg_list = verticall.pairwise.get_arg_list(Args(part='1/1'), assemblies, reference)
    assert len(arg_list) == 8


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


def test_check_assemblies():
    assemblies = [('a', pathlib.Path('test/test_pairwise/assemblies/a.fasta'))]
    verticall.pairwise.check_assemblies(assemblies, None)

    assemblies = [('b', pathlib.Path('test/test_pairwise/assemblies/b.fasta.gz'))]
    with pytest.raises(SystemExit) as e:
        verticall.pairwise.check_assemblies(assemblies, None)
    assert 'duplicate contig names' in str(e.value)

    assemblies = [('c', pathlib.Path('test/test_pairwise/assemblies/c.fna'))]
    with pytest.raises(SystemExit) as e:
        verticall.pairwise.check_assemblies(assemblies, None)
    assert 'ambiguous' in str(e.value)

