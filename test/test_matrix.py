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

import pytest

import verticall.matrix


def test_welcome_message(capsys):
    verticall.matrix.welcome_message()
    _, err = capsys.readouterr()
    assert 'Verticall matrix' in err


def test_finished_message(capsys):
    verticall.matrix.finished_message()
    _, err = capsys.readouterr()
    assert 'Finished' in err


def test_jukes_cantor():
    assert verticall.matrix.jukes_cantor(0.0) == 0.0
    assert verticall.matrix.jukes_cantor(0.9) == 25.0
    assert verticall.matrix.jukes_cantor(None) is None


def test_jukes_cantor_correction():
    sample_names = ['a', 'b']
    distances = {('a', 'a'): 0.0, ('a', 'b'): 0.2,
                 ('b', 'a'): 0.1, ('b', 'b'): 0.0}
    verticall.matrix.jukes_cantor_correction(distances, sample_names)
    assert distances[('a', 'a')] == pytest.approx(0.0)
    assert distances[('a', 'b')] == pytest.approx(0.23261619622788)
    assert distances[('b', 'a')] == pytest.approx(0.107325632730505)
    assert distances[('b', 'b')] == pytest.approx(0.0)


def test_make_symmetrical_1():
    sample_names = ['a', 'b']
    distances = {('a', 'a'): 0.0, ('a', 'b'): 0.2,
                 ('b', 'a'): 0.1, ('b', 'b'): 0.0}
    verticall.matrix.make_symmetrical(distances, sample_names)
    assert distances[('a', 'a')] == pytest.approx(0.0)
    assert distances[('a', 'b')] == pytest.approx(0.15)
    assert distances[('b', 'a')] == pytest.approx(0.15)
    assert distances[('b', 'b')] == pytest.approx(0.0)


def test_make_symmetrical_2():
    sample_names = ['a', 'b']
    distances = {('a', 'a'): 0.0, ('a', 'b'): None,
                 ('b', 'a'): 0.1, ('b', 'b'): 0.0}
    verticall.matrix.make_symmetrical(distances, sample_names)
    assert distances[('a', 'a')] == pytest.approx(0.0)
    assert distances[('a', 'b')] == pytest.approx(0.1)
    assert distances[('b', 'a')] == pytest.approx(0.1)
    assert distances[('b', 'b')] == pytest.approx(0.0)


def test_make_symmetrical_3():
    sample_names = ['a', 'b']
    distances = {('a', 'a'): 0.0,  ('a', 'b'): 0.2,
                 ('b', 'a'): None, ('b', 'b'): 0.0}
    verticall.matrix.make_symmetrical(distances, sample_names)
    assert distances[('a', 'a')] == pytest.approx(0.0)
    assert distances[('a', 'b')] == pytest.approx(0.2)
    assert distances[('b', 'a')] == pytest.approx(0.2)
    assert distances[('b', 'b')] == pytest.approx(0.0)


def test_make_symmetrical_4():
    sample_names = ['a', 'b']
    distances = {('a', 'a'): 0.0,  ('a', 'b'): None,
                 ('b', 'a'): None, ('b', 'b'): 0.0}
    verticall.matrix.make_symmetrical(distances, sample_names)
    assert distances[('a', 'a')] == pytest.approx(0.0)
    assert distances[('a', 'b')] is None
    assert distances[('b', 'a')] is None
    assert distances[('b', 'b')] == pytest.approx(0.0)


def test_check_for_missing_distances():
    sample_names = ['a', 'b', 'c']
    distances = {('a', 'a'): 0.0,
                 ('b', 'a'): 0.1, ('b', 'b'): 0.0, ('b', 'c'): 0.2,
                                  ('c', 'b'): 0.1, ('c', 'c'): 0.0}
    verticall.matrix.check_for_missing_distances(distances, sample_names)
    assert distances[('a', 'a')] == 0.0
    assert distances[('a', 'b')] is None
    assert distances[('a', 'c')] is None
    assert distances[('b', 'a')] == 0.1
    assert distances[('b', 'b')] == 0.0
    assert distances[('b', 'c')] == 0.2
    assert distances[('c', 'a')] is None
    assert distances[('c', 'b')] == 0.1
    assert distances[('c', 'c')] == 0.0


def test_filter_names():
    all_names = ['a', 'b', 'c', 'd', 'e', 'f']
    assert verticall.matrix.filter_names(all_names, 'b,c') == ['b', 'c']
    assert verticall.matrix.filter_names(all_names, 'f,c') == ['c', 'f']
    assert verticall.matrix.filter_names(all_names, 'e') == ['e']
    assert verticall.matrix.filter_names(all_names, 'a,b,c,d,e,f') == ['a', 'b', 'c', 'd', 'e', 'f']
    with pytest.raises(SystemExit) as e:
        verticall.matrix.filter_names(all_names, 'a,b,c,q,e,f')
    assert 'could not find sample' in str(e.value)


def test_multi_distance():
    assert verticall.matrix.multi_distance(0.4, 0.3, 'first') == 0.4
    assert verticall.matrix.multi_distance(0.4, 0.5, 'first') == 0.4
    assert verticall.matrix.multi_distance(0.4, 0.3, 'low') == 0.3
    assert verticall.matrix.multi_distance(0.4, 0.5, 'low') == 0.4
    assert verticall.matrix.multi_distance(0.4, 0.3, 'high') == 0.4
    assert verticall.matrix.multi_distance(0.4, 0.5, 'high') == 0.5
    with pytest.raises(AssertionError):
        verticall.matrix.multi_distance(0.4, 0.5, 'bad')


def test_get_distance_from_line_parts():
    parts = ['a', 'b', '0.0002', '', '0.0001', 'not_a_num']
    assert verticall.matrix.get_distance_from_line_parts(parts, 2) == pytest.approx(0.0002)
    assert verticall.matrix.get_distance_from_line_parts(parts, 3) is None
    assert verticall.matrix.get_distance_from_line_parts(parts, 4) == pytest.approx(0.0001)
    with pytest.raises(SystemExit) as e:
        verticall.matrix.get_distance_from_line_parts(parts, 5)
    assert 'could not convert' in str(e.value)
    with pytest.raises(SystemExit) as e:
        verticall.matrix.get_distance_from_line_parts(parts, 6)
    assert 'column' in str(e.value)
