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

import verticall.tsv


def test_get_column_index_1():
    header_parts = ['assembly_a', 'assembly_b', 'alignment_count', 'aligned_fraction',
                    'window_size', 'window_count', 'mean_distance', 'median_distance',
                    'mass_peaks', 'peak_distance', 'mean_vertical_distance',
                    'median_vertical_distance']
    assert verticall.tsv.get_column_index(header_parts, 'mean_distance', 'filename') == 6
    assert verticall.tsv.get_column_index(header_parts, 'median_distance', 'filename') == 7
    assert verticall.tsv.get_column_index(header_parts, 'peak_distance', 'filename') == 9
    assert verticall.tsv.get_column_index(header_parts, 'mean_vertical_distance', 'filename') == 10
    assert verticall.tsv.get_column_index(header_parts,
                                          'median_vertical_distance','filename') == 11
    with pytest.raises(SystemExit) as e:
        verticall.tsv.get_column_index(header_parts, 'bad', 'filename')
    assert 'no column named' in str(e.value)


def test_get_column_index_2():
    header_parts = ['bad_column_name', 'assembly_b']
    with pytest.raises(SystemExit) as e:
        verticall.tsv.get_column_index(header_parts, 'mean_distance', 'filename')
    assert 'is not labelled' in str(e.value)
    header_parts = ['assembly_a', 'bad_column_name']
    with pytest.raises(SystemExit) as e:
        verticall.tsv.get_column_index(header_parts, 'mean_distance', 'filename')
    assert 'is not labelled' in str(e.value)


def test_split_region_str():
    assert verticall.tsv.split_region_str('contig_1:5-10') == ('contig_1', 5, 10)
    assert verticall.tsv.split_region_str('contig_2:123-654') == ('contig_2', 123, 654)

    with pytest.raises(SystemExit) as e:
        verticall.tsv.split_region_str('contig_2:abc-def')
    assert 'not correctly formatted' in str(e.value)


def test_get_start_end():
    assert verticall.tsv.get_start_end('contig_1:5-10') == (5, 10)
    assert verticall.tsv.get_start_end('contig_2:123-654') == (123, 654)

    with pytest.raises(SystemExit) as e:
        verticall.tsv.get_start_end('contig_2:abc-def')
    assert 'not correctly formatted' in str(e.value)


def test_check_header_for_assembly_a_regions():
    header_parts = ['assembly_a', 'assembly_b', 'other', 'stuff',
                    'assembly_a_vertical_regions', 'assembly_a_horizontal_regions',
                    'assembly_a_unaligned_regions']
    verticall.tsv.check_header_for_assembly_a_regions(header_parts, 'filename')

    header_parts = ['assembly_a', 'assembly_b', 'other', 'stuff',
                    'assembly_a_horizontal_regions', 'assembly_a_unaligned_regions']
    with pytest.raises(SystemExit) as e:
        verticall.tsv.check_header_for_assembly_a_regions(header_parts, 'filename')
    assert 'no column named' in str(e.value)

    header_parts = ['assembly_a', 'assembly_b', 'other', 'stuff',
                    'assembly_a_vertical_regions', 'assembly_a_unaligned_regions']
    with pytest.raises(SystemExit) as e:
        verticall.tsv.check_header_for_assembly_a_regions(header_parts, 'filename')
    assert 'no column named' in str(e.value)

    header_parts = ['assembly_a', 'assembly_b', 'other', 'stuff',
                    'assembly_a_vertical_regions', 'assembly_a_horizontal_regions']
    with pytest.raises(SystemExit) as e:
        verticall.tsv.check_header_for_assembly_a_regions(header_parts, 'filename')
    assert 'no column named' in str(e.value)
