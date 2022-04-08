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


def test_jukes_cantor_correction():
    sample_names = ['a', 'b']
    distances = {('a', 'a'): 0.0, ('a', 'b'): 0.2,
                 ('b', 'a'): 0.1, ('b', 'b'): 0.0}
    verticall.matrix.jukes_cantor_correction(distances, sample_names)
    assert distances[('a', 'a')] == pytest.approx(0.0)
    assert distances[('a', 'b')] == pytest.approx(0.23261619622788)
    assert distances[('b', 'a')] == pytest.approx(0.107325632730505)
    assert distances[('b', 'b')] == pytest.approx(0.0)


def test_make_symmetrical():
    sample_names = ['a', 'b']
    distances = {('a', 'a'): 0.0, ('a', 'b'): 0.2,
                 ('b', 'a'): 0.1, ('b', 'b'): 0.0}
    verticall.matrix.make_symmetrical(distances, sample_names)
    assert distances[('a', 'a')] == pytest.approx(0.0)
    assert distances[('a', 'b')] == pytest.approx(0.15)
    assert distances[('b', 'a')] == pytest.approx(0.15)
    assert distances[('b', 'b')] == pytest.approx(0.0)


def test_get_column_index():
    header_parts = ['assembly_a', 'assembly_b', 'alignment_count', 'aligned_fraction',
                    'window_size', 'window_count', 'mean_distance', 'median_distance',
                    'mass_peaks', 'peak_distance', 'mean_vertical_distance',
                    'median_vertical_distance']
    assert verticall.matrix.get_column_index(header_parts, 'mean', 'filename') == 6
    assert verticall.matrix.get_column_index(header_parts, 'median', 'filename') == 7
    assert verticall.matrix.get_column_index(header_parts, 'peak', 'filename') == 9
    assert verticall.matrix.get_column_index(header_parts, 'mean_vertical', 'filename') == 10
    assert verticall.matrix.get_column_index(header_parts, 'median_vertical', 'filename') == 11


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
