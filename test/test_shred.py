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

import pathlib

import phylo.shred


def test_trim_seq():
    assert phylo.shred.trim_seq('ACGTACGTACGTACGT', 16) == 'ACGTACGTACGTACGT'
    assert phylo.shred.trim_seq('ACGTACGTACGTACGT', 15) == 'ACGTACGTACGTACG'
    assert phylo.shred.trim_seq('ACGTACGTACGTACGT', 14) ==  'CGTACGTACGTACG'
    assert phylo.shred.trim_seq('ACGTACGTACGTACGT', 13) ==  'CGTACGTACGTAC'
    assert phylo.shred.trim_seq('ACGTACGTACGTACGT', 12) ==   'GTACGTACGTAC'
    assert phylo.shred.trim_seq('ACGTACGTACGTACGT', 11) ==   'GTACGTACGTA'
    assert phylo.shred.trim_seq('ACGTACGTACGTACGT', 10) ==    'TACGTACGTA'
    assert phylo.shred.trim_seq('ACGTACGTACGTACGT', 9) ==     'TACGTACGT'
    assert phylo.shred.trim_seq('ACGTACGTACGTACGT', 8) ==      'ACGTACGT'
    assert phylo.shred.trim_seq('ACGTACGTACGTACGT', 7) ==      'ACGTACG'
    assert phylo.shred.trim_seq('ACGTACGTACGTACGT', 6) ==       'CGTACG'
    assert phylo.shred.trim_seq('ACGTACGTACGTACGT', 5) ==       'CGTAC'
    assert phylo.shred.trim_seq('ACGTACGTACGTACGT', 4) ==        'GTAC'
    assert phylo.shred.trim_seq('ACGTACGTACGTACGT', 3) ==        'GTA'
    assert phylo.shred.trim_seq('ACGTACGTACGTACGT', 2) ==         'TA'
    assert phylo.shred.trim_seq('ACGTACGTACGTACGT', 1) ==         'T'
    assert phylo.shred.trim_seq('ACGTACGTACGTACGT', 0) ==          ''


def test_shred_sequence():
    seq = 'CGACGCAGATCGACGCTAGC'
    assert phylo.shred.shred_sequence(seq, 10, 0) == ['CGACGCAGAT', 'CGACGCTAGC']
    assert phylo.shred.shred_sequence(seq, 5, 0) == ['CGACG', 'CAGAT', 'CGACG', 'CTAGC']
    assert phylo.shred.shred_sequence(seq, 10, 5) == ['CGACGCAGAT', 'CAGATCGACG', 'CGACGCTAGC']
    assert phylo.shred.shred_sequence(seq, 10, 5) == ['CGACGCAGAT', 'CAGATCGACG', 'CGACGCTAGC']
    assert phylo.shred.shred_sequence(seq, 10, 9) == ['CGACGCAGAT', 'GACGCAGATC', 'ACGCAGATCG',
                                                      'CGCAGATCGA', 'GCAGATCGAC', 'CAGATCGACG',
                                                      'AGATCGACGC', 'GATCGACGCT', 'ATCGACGCTA',
                                                      'TCGACGCTAG', 'CGACGCTAGC']
    assert phylo.shred.shred_sequence(seq, 16, 0) == ['ACGCAGATCGACGCTA']
    assert phylo.shred.shred_sequence(seq, 16, 11) == ['ACGCAGATCGACGCTA']
    assert phylo.shred.shred_sequence(seq, 16, 12) == ['CGACGCAGATCGACGC', 'GCAGATCGACGCTAGC']
    assert phylo.shred.shred_sequence(seq, 5, 3) == ['CGACG', 'ACGCA', 'GCAGA', 'AGATC', 'ATCGA',
                                                     'CGACG', 'ACGCT', 'GCTAG']
    assert phylo.shred.shred_sequence(seq, 6, 1) == ['ACGCAG', 'GATCGA', 'ACGCTA']


def test_get_name_from_filename():
    assert phylo.shred.get_name_from_filename(pathlib.Path('abc/123/name.fasta')) == 'name'
    assert phylo.shred.get_name_from_filename(pathlib.Path('abc/123/name.fasta.gz')) == 'name'
    assert phylo.shred.get_name_from_filename(pathlib.Path('abc/123/name.gfa')) == 'name'
    assert phylo.shred.get_name_from_filename(pathlib.Path('abc/123/name.gfa.gz')) == 'name'
