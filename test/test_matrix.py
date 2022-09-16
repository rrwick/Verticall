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


def test_jukes_cantor_correction_1():
    sample_names = ['a', 'b']
    distances = {('a', 'a'): 0.0, ('a', 'b'): 0.2,
                 ('b', 'a'): 0.1, ('b', 'b'): 0.0}
    verticall.matrix.jukes_cantor_correction(distances, sample_names)
    assert distances[('a', 'a')] == pytest.approx(0.0)
    assert distances[('a', 'b')] == pytest.approx(0.23261619622788)
    assert distances[('b', 'a')] == pytest.approx(0.107325632730505)
    assert distances[('b', 'b')] == pytest.approx(0.0)


def test_jukes_cantor_correction_2():
    sample_names = ['a', 'b']
    distances = {('a', 'a'): 0.0, ('a', 'b'): 0.2,
                 ('b', 'a'): 0.1}
    verticall.matrix.jukes_cantor_correction(distances, sample_names)
    assert distances[('a', 'a')] == pytest.approx(0.0)
    assert distances[('a', 'b')] == pytest.approx(0.23261619622788)
    assert distances[('b', 'a')] == pytest.approx(0.107325632730505)


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


def test_make_symmetrical_5():
    sample_names = ['a', 'b']
    distances = {('a', 'a'): 0.0,  ('a', 'b'): 0.2,
                 ('b', 'b'): 0.0}
    verticall.matrix.make_symmetrical(distances, sample_names)
    assert distances[('a', 'a')] == pytest.approx(0.0)
    assert distances[('a', 'b')] == pytest.approx(0.2)
    assert distances[('b', 'a')] == pytest.approx(0.2)
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


def test_include_names():
    all_names = ['a', 'b', 'c', 'd', 'e', 'f']
    assert verticall.matrix.include_names(all_names, 'b,c') == ['b', 'c']
    assert verticall.matrix.include_names(all_names, 'f,c') == ['c', 'f']
    assert verticall.matrix.include_names(all_names, 'e') == ['e']
    assert verticall.matrix.include_names(all_names, 'a,b,c,d,e,f') == ['a', 'b', 'c', 'd', 'e', 'f']
    with pytest.raises(SystemExit) as e:
        verticall.matrix.include_names(all_names, 'a,b,c,q,e,f')
    assert 'could not find sample' in str(e.value)


def test_exclude_names():
    all_names = ['a', 'b', 'c', 'd', 'e', 'f']
    assert verticall.matrix.exclude_names(all_names, 'b,c') == ['a', 'd', 'e', 'f']
    assert verticall.matrix.exclude_names(all_names, 'f,c') == ['a', 'b', 'd', 'e']
    assert verticall.matrix.exclude_names(all_names, 'e') == ['a', 'b', 'c', 'd', 'f']
    with pytest.raises(SystemExit) as e:
        verticall.matrix.exclude_names(all_names, 'a,b,c,d,e,f')
    assert 'all samples have been excluded' in str(e.value)
    with pytest.raises(SystemExit) as e:
        verticall.matrix.exclude_names(all_names, 'a,b,c,q')
    assert 'could not find sample' in str(e.value)


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


def get_truth_mean_distances():
    return {('INF001', 'INF001'): [0.0], ('INF002', 'INF002'): [0.0], ('INF003', 'INF003'): [0.0],
            ('INF004', 'INF004'): [0.0], ('INF005', 'INF005'): [0.0], ('INF013', 'INF013'): [0.0],
            ('INF062', 'INF062'): [0.0], ('INF097', 'INF097'): [0.0],
            ('INF001', 'INF002'): [0.053168105], ('INF001', 'INF003'): [0.052912422],
            ('INF001', 'INF004'): [0.053534657], ('INF001', 'INF005'): [0.008425837],
            ('INF001', 'INF013'): [0.052854770], ('INF001', 'INF062'): [0.052212158],
            ('INF001', 'INF097'): [0.052881955], ('INF002', 'INF001'): [0.053200876],
            ('INF002', 'INF003'): [0.007447296], ('INF002', 'INF004'): [0.007081810],
            ('INF002', 'INF005'): [0.052245664], ('INF002', 'INF013'): [0.007337412],
            ('INF002', 'INF062'): [0.008255233], ('INF002', 'INF097'): [0.007536904],
            ('INF003', 'INF001'): [0.052891430], ('INF003', 'INF002'): [0.007537941],
            ('INF003', 'INF004'): [0.008005609], ('INF003', 'INF005'): [0.053103535],
            ('INF003', 'INF013'): [0.007065224], ('INF003', 'INF062'): [0.007853800],
            ('INF003', 'INF097'): [0.007284896], ('INF004', 'INF001'): [0.053639524],
            ('INF004', 'INF002'): [0.007110545], ('INF004', 'INF003'): [0.007961059],
            ('INF004', 'INF005'): [0.053485734], ('INF004', 'INF013'): [0.007115893],
            ('INF004', 'INF062'): [0.007946477], ('INF004', 'INF097'): [0.007379517],
            ('INF005', 'INF001'): [0.008435388], ('INF005', 'INF002'): [0.052337464],
            ('INF005', 'INF003'): [0.053162877], ('INF005', 'INF004'): [0.053480168],
            ('INF005', 'INF013'): [0.052818698], ('INF005', 'INF062'): [0.052736855],
            ('INF005', 'INF097'): [0.052507758], ('INF013', 'INF001'): [0.052804671],
            ('INF013', 'INF002'): [0.007396860], ('INF013', 'INF003'): [0.007089594],
            ('INF013', 'INF004'): [0.007106779], ('INF013', 'INF005'): [0.052822142],
            ('INF013', 'INF062'): [0.003480056, 0.003480056],
            ('INF013', 'INF097'): [0.003483088, 0.003483088],
            ('INF062', 'INF001'): [0.052217258], ('INF062', 'INF002'): [0.008228569],
            ('INF062', 'INF003'): [0.007815157], ('INF062', 'INF004'): [0.007936269],
            ('INF062', 'INF005'): [0.052756774],
            ('INF062', 'INF013'): [0.003496058, 0.003496058],
            ('INF062', 'INF097'): [0.003717277, 0.003717277],
            ('INF097', 'INF001'): [0.052848589], ('INF097', 'INF002'): [0.007527977],
            ('INF097', 'INF003'): [0.007304839], ('INF097', 'INF004'): [0.007413049],
            ('INF097', 'INF005'): [0.052586871],
            ('INF097', 'INF013'): [0.003470136, 0.003470136],
            ('INF097', 'INF062'): [0.003736065, 0.003736065]}


def get_truth_median_vertical_window_distances():
    return {('INF001', 'INF001'): [0.0], ('INF002', 'INF002'): [0.0], ('INF003', 'INF003'): [0.0],
            ('INF004', 'INF004'): [0.0], ('INF005', 'INF005'): [0.0], ('INF013', 'INF013'): [0.0],
            ('INF062', 'INF062'): [0.0], ('INF097', 'INF097'): [0.0],
            ('INF001', 'INF002'): [0.050198172], ('INF001', 'INF003'): [0.050200211],
            ('INF001', 'INF004'): [0.050569207], ('INF001', 'INF005'): [0.006463670],
            ('INF001', 'INF013'): [0.050372656], ('INF001', 'INF062'): [0.049926028],
            ('INF001', 'INF097'): [0.050067421], ('INF002', 'INF001'): [0.050141234],
            ('INF002', 'INF003'): [0.005623070], ('INF002', 'INF004'): [0.005556619],
            ('INF002', 'INF005'): [0.049044106], ('INF002', 'INF013'): [0.005438631],
            ('INF002', 'INF062'): [0.005596103], ('INF002', 'INF097'): [0.005513849],
            ('INF003', 'INF001'): [0.050165005], ('INF003', 'INF002'): [0.005621025],
            ('INF003', 'INF004'): [0.005610997], ('INF003', 'INF005'): [0.050056432],
            ('INF003', 'INF013'): [0.005220400], ('INF003', 'INF062'): [0.005387198],
            ('INF003', 'INF097'): [0.005375850], ('INF004', 'INF001'): [0.050532858],
            ('INF004', 'INF002'): [0.005557945], ('INF004', 'INF003'): [0.005618168],
            ('INF004', 'INF005'): [0.049781057], ('INF004', 'INF013'): [0.005449283],
            ('INF004', 'INF062'): [0.005514706], ('INF004', 'INF097'): [0.005517577],
            ('INF005', 'INF001'): [0.006463747], ('INF005', 'INF002'): [0.049095533],
            ('INF005', 'INF003'): [0.050026761], ('INF005', 'INF004'): [0.049735624],
            ('INF005', 'INF013'): [0.049800733], ('INF005', 'INF062'): [0.049759494],
            ('INF005', 'INF097'): [0.049547786], ('INF013', 'INF001'): [0.050320362],
            ('INF013', 'INF002'): [0.005445456], ('INF013', 'INF003'): [0.005220464],
            ('INF013', 'INF004'): [0.005447387], ('INF013', 'INF005'): [0.049834326],
            ('INF013', 'INF062'): [0.000093611, 0.005251469],
            ('INF013', 'INF097'): [0.000090086, 0.004955741],
            ('INF062', 'INF001'): [0.049819225], ('INF062', 'INF002'): [0.005599151],
            ('INF062', 'INF003'): [0.005386446], ('INF062', 'INF004'): [0.005516905],
            ('INF062', 'INF005'): [0.049709677],
            ('INF062', 'INF013'): [0.000093620, 0.005251208],
            ('INF062', 'INF097'): [0.000098037, 0.005210406],
            ('INF097', 'INF001'): [0.050000893], ('INF097', 'INF002'): [0.005519458],
            ('INF097', 'INF003'): [0.005382644], ('INF097', 'INF004'): [0.005517577],
            ('INF097', 'INF005'): [0.049559676],
            ('INF097', 'INF013'): [0.000090078, 0.004955854],
            ('INF097', 'INF062'): [0.000098039, 0.005209505]}


def test_load_tsv_file_1():
    in_tsv = pathlib.Path('test/test_matrix/pairwise.tsv')
    distances, sample_names = verticall.matrix.load_tsv_file(in_tsv, 'mean')
    truth = get_truth_mean_distances()
    for a in sample_names:
        for b in sample_names:
            assert distances[(a, b)] == pytest.approx(truth[(a, b)])


def test_load_tsv_file_2():
    in_tsv = pathlib.Path('test/test_matrix/pairwise.tsv')
    distances, sample_names = verticall.matrix.load_tsv_file(in_tsv, 'median_vertical_window')
    truth = get_truth_median_vertical_window_distances()
    for a in sample_names:
        for b in sample_names:
            assert distances[(a, b)] == pytest.approx(truth[(a, b)])


def test_load_tsv_file_3():
    in_tsv = pathlib.Path('test/test_matrix/pairwise.tsv.gz')
    distances, sample_names = verticall.matrix.load_tsv_file(in_tsv, 'mean')
    truth = get_truth_mean_distances()
    for a in sample_names:
        for b in sample_names:
            assert distances[(a, b)] == pytest.approx(truth[(a, b)])


def test_load_tsv_file_4():
    in_tsv = pathlib.Path('test/test_matrix/pairwise.tsv.gz')
    distances, sample_names = verticall.matrix.load_tsv_file(in_tsv, 'median_vertical_window')
    truth = get_truth_median_vertical_window_distances()
    for a in sample_names:
        for b in sample_names:
            assert distances[(a, b)] == pytest.approx(truth[(a, b)])


def test_resolve_multi_distances_1():
    in_tsv = pathlib.Path('test/test_matrix/pairwise.tsv.gz')
    distances, sample_names = verticall.matrix.load_tsv_file(in_tsv, 'median_vertical_window')
    distances, sample_names = verticall.matrix.resolve_multi_distances(distances, sample_names,
                                                                       'first')
    assert sample_names == ['INF001', 'INF002', 'INF003', 'INF004',
                            'INF005', 'INF013', 'INF062', 'INF097']
    truth = get_truth_median_vertical_window_distances()
    for a in sample_names:
        for b in sample_names:
            if a == b:
                assert distances[(a, b)] == pytest.approx(truth[(a, b)][0])


def test_resolve_multi_distances_2():
    in_tsv = pathlib.Path('test/test_matrix/pairwise.tsv.gz')
    distances, sample_names = verticall.matrix.load_tsv_file(in_tsv, 'median_vertical_window')
    distances, sample_names = verticall.matrix.resolve_multi_distances(distances, sample_names,
                                                                       'low')
    assert sample_names == ['INF001', 'INF002', 'INF003', 'INF004',
                            'INF005', 'INF013', 'INF062', 'INF097']
    truth = get_truth_median_vertical_window_distances()
    for a in sample_names:
        for b in sample_names:
            assert distances[(a, b)] == pytest.approx(min(truth[(a, b)]))


def test_resolve_multi_distances_3():
    in_tsv = pathlib.Path('test/test_matrix/pairwise.tsv.gz')
    distances, sample_names = verticall.matrix.load_tsv_file(in_tsv, 'median_vertical_window')
    distances, sample_names = verticall.matrix.resolve_multi_distances(distances, sample_names,
                                                                       'high')
    assert sample_names == ['INF001', 'INF002', 'INF003', 'INF004',
                            'INF005', 'INF013', 'INF062', 'INF097']
    truth = get_truth_median_vertical_window_distances()
    for a in sample_names:
        for b in sample_names:
            assert distances[(a, b)] == pytest.approx(max(truth[(a, b)]))


def test_resolve_multi_distances_4():
    in_tsv = pathlib.Path('test/test_matrix/pairwise.tsv.gz')
    distances, sample_names = verticall.matrix.load_tsv_file(in_tsv, 'median_vertical_window')
    distances, sample_names = verticall.matrix.resolve_multi_distances(distances, sample_names,
                                                                       'exclude')
    assert sample_names == ['INF001', 'INF002', 'INF003', 'INF004', 'INF005']
    truth = get_truth_median_vertical_window_distances()
    for a in sample_names:
        for b in sample_names:
            assert len(truth[(a, b)]) == 1
            assert distances[(a, b)] == pytest.approx(truth[(a, b)][0])


def test_resolve_multi_distances_5():
    in_tsv = pathlib.Path('test/test_matrix/pairwise.tsv.gz')
    distances, sample_names = verticall.matrix.load_tsv_file(in_tsv, 'median_vertical_window')
    sample_names = ['INF001', 'INF002', 'INF003', 'INF004', 'INF005']
    del distances[('INF013', 'INF062')]
    del distances[('INF013', 'INF097')]
    del distances[('INF062', 'INF013')]
    del distances[('INF062', 'INF097')]
    del distances[('INF097', 'INF013')]
    del distances[('INF097', 'INF062')]
    distances, sample_names = verticall.matrix.resolve_multi_distances(distances, sample_names,
                                                                       'first')
    truth = get_truth_median_vertical_window_distances()
    for a in sample_names:
        for b in sample_names:
            assert len(truth[(a, b)]) == 1
            assert distances[(a, b)] == pytest.approx(truth[(a, b)][0])


def test_resolve_multi_distances_6():
    in_tsv = pathlib.Path('test/test_matrix/pairwise.tsv.gz')
    distances, sample_names = verticall.matrix.load_tsv_file(in_tsv, 'median_vertical_window')
    with pytest.raises(AssertionError):
        verticall.matrix.resolve_multi_distances(distances, sample_names, 'not_an_option')


def test_save_matrix_1():
    sample_names = ['a', 'b', 'c']
    distances = {('a', 'a'): 0.0, ('a', 'b'): 0.1, ('a', 'c'): 0.2,
                 ('b', 'a'): 0.1, ('b', 'b'): 0.0, ('b', 'c'): 0.3,
                 ('c', 'a'): 0.2, ('c', 'b'): 0.3, ('c', 'c'): 0.0}
    with tempfile.TemporaryDirectory() as temp_dir:
        matrix_filename = pathlib.Path(temp_dir) / 'matrix.phylip'
        verticall.matrix.save_matrix(matrix_filename, distances, sample_names)
        with open(matrix_filename, 'rt') as f:
            matrix = f.read()
            assert matrix == ('3\n'
                              'a\t0.000000000\t0.100000000\t0.200000000\n'
                              'b\t0.100000000\t0.000000000\t0.300000000\n'
                              'c\t0.200000000\t0.300000000\t0.000000000\n')


def test_save_matrix_2():
    sample_names = ['a', 'b', 'c']
    distances = {('a', 'a'): 0.0, ('a', 'b'): 0.1, ('a', 'c'): 0.2,
                 ('b', 'a'): 0.1, ('b', 'b'): 0.0, ('b', 'c'): 0.3,
                 ('c', 'a'): 0.2, ('c', 'b'): 0.3, ('c', 'c'): 0.0}
    with tempfile.TemporaryDirectory() as temp_dir:
        matrix_filename = pathlib.Path(temp_dir) / 'matrix.phylip'
        verticall.matrix.save_matrix(matrix_filename, distances, sample_names, silent=True)
        with open(matrix_filename, 'rt') as f:
            matrix = f.read()
            assert matrix == ('3\n'
                              'a\t0.000000000\t0.100000000\t0.200000000\n'
                              'b\t0.100000000\t0.000000000\t0.300000000\n'
                              'c\t0.200000000\t0.300000000\t0.000000000\n')


def test_save_matrix_3():
    sample_names = ['a', 'b', 'c']
    distances = {('a', 'a'): 0.0, ('a', 'b'): 0.1, ('a', 'c'): 0.2,
                 ('b', 'a'): 0.1, ('b', 'b'): 0.0,
                 ('c', 'a'): None, ('c', 'b'): 0.3, ('c', 'c'): 0.0}
    with tempfile.TemporaryDirectory() as temp_dir:
        matrix_filename = pathlib.Path(temp_dir) / 'matrix.phylip'
        verticall.matrix.save_matrix(matrix_filename, distances, sample_names)
        with open(matrix_filename, 'rt') as f:
            matrix = f.read()
            assert matrix == ('3\n'
                              'a\t0.000000000\t0.100000000\t0.200000000\n'
                              'b\t0.100000000\t0.000000000\t\n'
                              'c\t\t0.300000000\t0.000000000\n')
