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
    sequences = {'A': 'AGCTACGACCTA',
                 'B': 'AGCTACG-CCTA',
                 'C': 'AGCNNNGACCTA'}
    assert verticall.mask.drop_empty_positions(sequences) == {'A': 'AGCTACGACCTA',
                                                              'B': 'AGCTACG-CCTA',
                                                              'C': 'AGCNNNGACCTA'}
    sequences = {'A': 'AGCTACG-CCTA',
                 'B': 'AGCTACG-CCTA',
                 'C': 'AGCNNNG-CCTA'}
    assert verticall.mask.drop_empty_positions(sequences) == {'A': 'AGCTACGCCTA',
                                                              'B': 'AGCTACGCCTA',
                                                              'C': 'AGCNNNGCCTA'}
    sequences = {'A': 'AGNT-CNACN-A',
                 'B': 'AGNT-C--C--A',
                 'C': 'AGNNNNG-C--A'}
    assert verticall.mask.drop_empty_positions(sequences) == {'A': 'AGTCNACA',
                                                              'B': 'AGTC--CA',
                                                              'C': 'AGNNG-CA'}


def test_drop_invariant_positions():
    sequences = {'A': 'AGCTACGACCTA',
                 'B': 'AGCAACGACGTA',
                 'C': 'AGCTACGCCCTA'}
    assert verticall.mask.drop_invariant_positions(sequences) == {'A': 'TAC',
                                                                  'B': 'AAG',
                                                                  'C': 'TCC'}
    sequences = {'A': 'AGCTACGACCTA',
                 'B': 'TCGACTGACGAC',
                 'C': 'ACGACTACGACG'}
    assert verticall.mask.drop_invariant_positions(sequences) == {'A': 'AGCTACGACCTA',
                                                                  'B': 'TCGACTGACGAC',
                                                                  'C': 'ACGACTACGACG'}
    sequences = {'A': 'AGCTACGACCTA',
                 'B': 'AGCTACGACCTA',
                 'C': 'AGCTACGACCTA'}
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
