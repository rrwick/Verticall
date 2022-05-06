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

import verticall.mask


def test_welcome_message(capsys):
    Args = collections.namedtuple('Args', ['in_tsv', 'in_alignment', 'out_alignment'])

    verticall.mask.welcome_message(Args(in_tsv='in.tsv', in_alignment='in.fasta',
                                        out_alignment='out.fasta'))
    _, err = capsys.readouterr()
    assert 'Verticall mask' in err


def test_finished_message(capsys):
    verticall.mask.finished_message()
    _, err = capsys.readouterr()
    assert 'Finished' in err


def test_get_start_end():
    assert verticall.mask.get_start_end('a:1-2') == (1, 2)
    assert verticall.mask.get_start_end('b:10-20') == (10, 20)


def test_get_alignment_positions():
    #   aligned positions: 01234567
    #         aligned seq: ACGATCGA
    # unaligned positions: 01234567
    assert verticall.mask.get_alignment_positions('ACGATCGA', 8) == \
           {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8}

    #   aligned positions: 012345678
    #         aligned seq: ACGA-TCGA
    # unaligned positions: 0123 4567
    assert verticall.mask.get_alignment_positions('ACGA-TCGA', 8) == \
           {0: 0, 1: 1, 2: 2, 3: 3, 4: 5, 5: 6, 6: 7, 7: 8, 8: 9}

    #   aligned positions: 0123456789
    #         aligned seq: A-CGA-TCGA
    # unaligned positions: 0 123 4567
    assert verticall.mask.get_alignment_positions('A-CGA-TCGA', 8) == \
           {0: 0, 1: 2, 2: 3, 3: 4, 4: 6, 5: 7, 6: 8, 7: 9, 8: 10}

    with pytest.raises(SystemExit) as e:
        verticall.mask.get_alignment_positions('A--GA--CGA', 8)
    assert 'length of reference sequence' in str(e.value)


def test_drop_empty_positions():
    sequences = {'A': 'agcTACGACCTA',
                 'B': 'agcTACG-CcTA',
                 'C': 'agcNNNGACcTA'}
    assert verticall.mask.drop_empty_positions(sequences) == {'A': 'agcTACGACCTA',
                                                              'B': 'agcTACG-CcTA',
                                                              'C': 'agcNNNGACcTA'}
    sequences = {'A': 'AGCtACG-CCTA',
                 'B': 'AGCtACg-CCTA',
                 'C': 'AGCNNNG-CCTA'}
    assert verticall.mask.drop_empty_positions(sequences) == {'A': 'AGCtACGCCTA',
                                                              'B': 'AGCtACgCCTA',
                                                              'C': 'AGCNNNGCCTA'}
    sequences = {'A': 'AGNT-CNACN-A',
                 'B': 'AGNT-C--C--A',
                 'C': 'AGNNNNG-C--A'}
    assert verticall.mask.drop_empty_positions(sequences) == {'A': 'AGTCNACA',
                                                              'B': 'AGTC--CA',
                                                              'C': 'AGNNG-CA'}


def test_drop_invariant_positions():
    sequences = {'A': 'AGCTACGAcctA',
                 'B': 'aGCaACGACGtA',
                 'C': 'AGCTACGCcCtA'}
    assert verticall.mask.drop_invariant_positions(sequences) == {'A': 'TAc',
                                                                  'B': 'aAG',
                                                                  'C': 'TCC'}
    sequences = {'A': 'AGCTACGACCTA',
                 'B': 'TCGACTGACGAC',
                 'C': 'ACGACTACGACG'}
    assert verticall.mask.drop_invariant_positions(sequences) == {'A': 'AGCTACGACCTA',
                                                                  'B': 'TCGACTGACGAC',
                                                                  'C': 'ACGACTACGACG'}
    sequences = {'A': 'agCTACGACcta',
                 'B': 'agCTacgACCTA',
                 'C': 'agCTACGACCTA'}
    assert verticall.mask.drop_invariant_positions(sequences) == {'A': '',
                                                                  'B': '',
                                                                  'C': ''}


def test_count_real_bases():
    assert verticall.mask.count_real_bases({'N'}) == 0
    assert verticall.mask.count_real_bases({'A'}) == 1
    assert verticall.mask.count_real_bases({'A', 'T'}) == 2
    assert verticall.mask.count_real_bases({'A', 'N'}) == 1
    assert verticall.mask.count_real_bases({'-', 'N'}) == 0
    assert verticall.mask.count_real_bases({'A', 'N', 'C', '-', 'G', 'X', 'T'}) == 4


def test_load_regions():
    in_tsv = pathlib.Path('test/test_mask/pairwise.tsv')

    data, ref_name, ref_length, sample_names = verticall.mask.load_regions(in_tsv, 'ref', 'exclude')
    assert ref_name == 'ref'
    assert ref_length == 3000
    assert sample_names == ['1', '2', '3']
    assert data['1'] == ([(0, 1000), (2000, 3000)], [(1000, 2000)], [])
    assert data['2'] == ([(0, 1100), (2100, 3000)], [(1100, 2000)], [(2000, 2100)])
    assert data['3'] == ([(0, 1200), (2200, 3000)], [(1200, 2000)], [(2000, 2200)])

    data, ref_name, ref_length, sample_names = verticall.mask.load_regions(in_tsv, 'ref', 'first')
    assert sample_names == ['1', '2', '3', '4']
    assert data['4'] == ([(0, 1000)], [(1000, 3000)], [])

    data, ref_name, ref_length, sample_names = verticall.mask.load_regions(in_tsv, 'ref', 'low')
    assert data['4'] == ([(1000, 2000)], [(0, 1000), (2000, 3000)], [])

    data, ref_name, ref_length, sample_names = verticall.mask.load_regions(in_tsv, 'ref', 'high')
    assert data['4'] == ([(2000, 3000)], [(0, 2000)], [])


def test_get_ref_length():
    data = {'1': ([(0, 1000), (2000, 5000)], [(1000, 2000)], []),
            '2': ([(0, 5000)], [], []),
            '3': ([], [], [(0, 5000)])}
    assert verticall.mask.get_ref_length(data) == 5000

    data = {'1': ([(0, 1000), (2000, 5000)], [(1000, 2000)], []),
            '2': ([(0, 6000)], [], []),
            '3': ([], [], [(0, 5000)])}
    with pytest.raises(SystemExit) as e:
        verticall.mask.get_ref_length(data)
    assert 'inconsistent' in str(e.value)
