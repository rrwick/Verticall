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


def test_finalise():
    sequences = {'A': 'AACtACG-CCTA',
                 'B': 'AGCtACg-CCGA',
                 'C': 'AGCNNNG-CCTA'}
    assert verticall.mask.finalise(sequences, False) == {'A': 'AACtACGCCTA',
                                                         'B': 'AGCtACgCCGA',
                                                         'C': 'AGCNNNGCCTA'}
    assert verticall.mask.finalise(sequences, True) == {'A': 'AT',
                                                        'B': 'GG',
                                                        'C': 'GT'}


def test_count_real_bases():
    assert verticall.mask.count_real_bases({'N'}) == 0
    assert verticall.mask.count_real_bases({'A'}) == 1
    assert verticall.mask.count_real_bases({'A', 'T'}) == 2
    assert verticall.mask.count_real_bases({'A', 'N'}) == 1
    assert verticall.mask.count_real_bases({'-', 'N'}) == 0
    assert verticall.mask.count_real_bases({'A', 'N', 'C', '-', 'G', 'X', 'T'}) == 4


def test_check_tsv_file_1():
    in_tsv = pathlib.Path('test/test_mask/pairwise.tsv')
    ref_name = verticall.mask.check_tsv_file(in_tsv, 'ref')
    assert ref_name == 'ref'


def test_check_tsv_file_2():
    # Fails to auto-determine reference name because there are two names in column 1.
    in_tsv = pathlib.Path('test/test_mask/pairwise.tsv')
    with pytest.raises(SystemExit) as e:
        verticall.mask.check_tsv_file(in_tsv, None)
    assert 'could not automatically determine the reference name' in str(e.value)


def test_check_tsv_file_3():
    # Succeeds in auto-determining reference name because there is only one name in column 1.
    in_tsv = pathlib.Path('test/test_mask/unambiguous_ref.tsv')
    ref_name = verticall.mask.check_tsv_file(in_tsv, None)
    assert ref_name == 'ref'


def test_load_regions_1():
    in_tsv = pathlib.Path('test/test_mask/pairwise.tsv')
    data, ref_name, ref_length, sample_names = verticall.mask.load_regions(in_tsv, 'ref', 'first')
    assert ref_name == 'ref'
    assert ref_length == 30
    assert sample_names == ['1', '2', '3', '4']
    assert data['1'] == ([(0, 10), (20, 30)], [(10, 20)], [])
    assert data['2'] == ([(0, 11), (21, 30)], [(11, 20)], [(20, 21)])
    assert data['3'] == ([(0, 12), (22, 30)], [(12, 20)], [(20, 22)])
    assert data['4'] == ([(0, 10)], [(10, 30)], [])


def test_load_regions_2():
    in_tsv = pathlib.Path('test/test_mask/pairwise.tsv')
    data, ref_name, ref_length, sample_names = verticall.mask.load_regions(in_tsv, 'ref', 'exclude')
    assert ref_name == 'ref'
    assert ref_length == 30
    assert sample_names == ['1', '2', '3']
    assert data['1'] == ([(0, 10), (20, 30)], [(10, 20)], [])
    assert data['2'] == ([(0, 11), (21, 30)], [(11, 20)], [(20, 21)])
    assert data['3'] == ([(0, 12), (22, 30)], [(12, 20)], [(20, 22)])
    assert '4' not in data


def test_load_regions_3():
    in_tsv = pathlib.Path('test/test_mask/pairwise.tsv')
    data, ref_name, ref_length, sample_names = verticall.mask.load_regions(in_tsv, 'ref', 'low')
    assert ref_name == 'ref'
    assert ref_length == 30
    assert data['1'] == ([(0, 10), (20, 30)], [(10, 20)], [])
    assert data['2'] == ([(0, 11), (21, 30)], [(11, 20)], [(20, 21)])
    assert data['3'] == ([(0, 12), (22, 30)], [(12, 20)], [(20, 22)])
    assert data['4'] == ([(10, 20)], [(0, 10), (20, 30)], [])


def test_load_regions_4():
    in_tsv = pathlib.Path('test/test_mask/pairwise.tsv')
    data, ref_name, ref_length, sample_names = verticall.mask.load_regions(in_tsv, 'ref', 'high')
    assert ref_name == 'ref'
    assert ref_length == 30
    assert data['1'] == ([(0, 10), (20, 30)], [(10, 20)], [])
    assert data['2'] == ([(0, 11), (21, 30)], [(11, 20)], [(20, 21)])
    assert data['3'] == ([(0, 12), (22, 30)], [(12, 20)], [(20, 22)])
    assert data['4'] == ([(20, 30)], [(0, 20)], [])


def test_load_regions_5():
    in_tsv = pathlib.Path('test/test_mask/pairwise.tsv')
    with pytest.raises(SystemExit) as e:
        verticall.mask.load_regions(in_tsv, 'bad_ref_name', 'first')
    assert 'reference-to-assembly pairwise comparisons found' in str(e.value)


def test_load_regions_6():
    in_tsv = pathlib.Path('test/test_mask/multicontig_ref.tsv')
    with pytest.raises(SystemExit) as e:
        verticall.mask.load_regions(in_tsv, 'ref', 'first')
    assert 'reference genome has more than one contig name' in str(e.value)


def test_load_pseudo_alignment_1():
    in_align = pathlib.Path('test/test_mask/alignment.fasta')
    sample_names = ['1', '2', '3', '4']
    sequences, sample_names = verticall.mask.load_pseudo_alignment(in_align, 'ref', sample_names)
    assert sample_names == ['1', '2', '3', '4']
    assert sequences['ref'] == 'GTACGCATCTCTTC--TCTGTAGCAATGAGAT'
    assert sequences['1'] ==   'GTAcnnatCTCTTCAATCTGTAGCAATGAGAT'
    assert sequences['2'] ==   'GTACGCATCTCTTC--TCtgtagcaATNNNAT'
    assert sequences['3'] ==   'GTACGCATctcttc--tcTGTAGCAATGA---'
    assert sequences['4'] ==   'GTA-----CTCTTC--TCTGTAgcaatgagAT'


def test_load_pseudo_alignment_2():
    in_align = pathlib.Path('test/test_mask/alignment.fasta')
    sample_names = ['A', 'B', 'C', 'D']
    with pytest.raises(SystemExit) as e:
        verticall.mask.load_pseudo_alignment(in_align, 'bad_ref_name', sample_names)
    assert 'could not find reference sequence' in str(e.value)


def test_load_pseudo_alignment_3():
    in_align = pathlib.Path('test/test_mask/empty.fasta')
    sample_names = ['A', 'B', 'C', 'D']
    with pytest.raises(SystemExit) as e:
        verticall.mask.load_pseudo_alignment(in_align, 'ref', sample_names)
    assert 'no sequences could be loaded' in str(e.value)


def test_load_pseudo_alignment_4():
    in_align = pathlib.Path('test/test_mask/alignment.fasta')
    sample_names = ['A', 'B', 'C', 'D']
    with pytest.raises(SystemExit) as e:
        verticall.mask.load_pseudo_alignment(in_align, 'ref', sample_names)
    assert 'no sample names in common' in str(e.value)


def test_load_pseudo_alignment_5():
    in_align = pathlib.Path('test/test_mask/different_lengths.fasta')
    sample_names = ['A', 'B', 'C', 'D']
    with pytest.raises(SystemExit) as e:
        verticall.mask.load_pseudo_alignment(in_align, 'ref', sample_names)
    assert 'must be the same length' in str(e.value)


def test_mask_sequences_1():
    in_tsv = pathlib.Path('test/test_mask/pairwise.tsv')
    in_align = pathlib.Path('test/test_mask/alignment.fasta')
    data, ref_name, ref_length, sample_names = verticall.mask.load_regions(in_tsv, 'ref', 'first')
    sequences, sample_names = verticall.mask.load_pseudo_alignment(in_align, ref_name, sample_names)
    masked_sequences = \
        verticall.mask.mask_sequences(data, sequences, ref_name, ref_length, sample_names,
                                      'N', '-', None, '#4859a0', '#c47e7e', '#c9c9c9')
    assert list(masked_sequences.keys()) == ['ref', '1', '2', '3', '4']
    assert masked_sequences['ref'] == 'GTACGCATCTCTTC--TCTGTAGCAATGAGAT'
    assert masked_sequences['1'] ==   'GTAcnnatCTNNNNNNNNNNNNGCAATGAGAT'
    assert masked_sequences['2'] ==   'GTACGCATCTCNNNNNNNNNNN-caATNNNAT'
    assert masked_sequences['3'] ==   'GTACGCATctctNNNNNNNNNN--AATGA---'
    assert masked_sequences['4'] ==   'GTA-----CTNNNNNNNNNNNNNNNNNNNNNN'


def test_mask_sequences_2():
    in_tsv = pathlib.Path('test/test_mask/pairwise.tsv')
    in_align = pathlib.Path('test/test_mask/alignment.fasta')
    data, ref_name, ref_length, sample_names = verticall.mask.load_regions(in_tsv, 'ref', 'exclude')
    sequences, sample_names = verticall.mask.load_pseudo_alignment(in_align, ref_name, sample_names)
    masked_sequences = \
        verticall.mask.mask_sequences(data, sequences, ref_name, ref_length, sample_names,
                                      'N', '-', None, '#4859a0', '#c47e7e', '#c9c9c9')
    assert list(masked_sequences.keys()) == ['ref', '1', '2', '3']
    assert masked_sequences['ref'] == 'GTACGCATCTCTTC--TCTGTAGCAATGAGAT'
    assert masked_sequences['1'] ==   'GTAcnnatCTNNNNNNNNNNNNGCAATGAGAT'
    assert masked_sequences['2'] ==   'GTACGCATCTCNNNNNNNNNNN-caATNNNAT'
    assert masked_sequences['3'] ==   'GTACGCATctctNNNNNNNNNN--AATGA---'


def test_mask_sequences_3():
    in_tsv = pathlib.Path('test/test_mask/pairwise.tsv')
    in_align = pathlib.Path('test/test_mask/alignment.fasta')
    data, ref_name, ref_length, sample_names = verticall.mask.load_regions(in_tsv, 'ref', 'low')
    sequences, sample_names = verticall.mask.load_pseudo_alignment(in_align, ref_name, sample_names)
    masked_sequences = \
        verticall.mask.mask_sequences(data, sequences, ref_name, ref_length, sample_names,
                                      'N', '-', None, '#4859a0', '#c47e7e', '#c9c9c9')
    assert list(masked_sequences.keys()) == ['ref', '1', '2', '3', '4']
    assert masked_sequences['ref'] == 'GTACGCATCTCTTC--TCTGTAGCAATGAGAT'
    assert masked_sequences['1'] ==   'GTAcnnatCTNNNNNNNNNNNNGCAATGAGAT'
    assert masked_sequences['2'] ==   'GTACGCATCTCNNNNNNNNNNN-caATNNNAT'
    assert masked_sequences['3'] ==   'GTACGCATctctNNNNNNNNNN--AATGA---'
    assert masked_sequences['4'] ==   'NNNNNNNNNNCTTC--TCTGTANNNNNNNNNN'


def test_mask_sequences_4():
    in_tsv = pathlib.Path('test/test_mask/pairwise.tsv')
    in_align = pathlib.Path('test/test_mask/alignment.fasta')
    data, ref_name, ref_length, sample_names = verticall.mask.load_regions(in_tsv, 'ref', 'high')
    sequences, sample_names = verticall.mask.load_pseudo_alignment(in_align, ref_name, sample_names)
    masked_sequences = \
        verticall.mask.mask_sequences(data, sequences, ref_name, ref_length, sample_names,
                                      'N', '-', None, '#4859a0', '#c47e7e', '#c9c9c9')
    assert list(masked_sequences.keys()) == ['ref', '1', '2', '3', '4']
    assert masked_sequences['ref'] == 'GTACGCATCTCTTC--TCTGTAGCAATGAGAT'
    assert masked_sequences['1'] ==   'GTAcnnatCTNNNNNNNNNNNNGCAATGAGAT'
    assert masked_sequences['2'] ==   'GTACGCATCTCNNNNNNNNNNN-caATNNNAT'
    assert masked_sequences['3'] ==   'GTACGCATctctNNNNNNNNNN--AATGA---'
    assert masked_sequences['4'] ==   'NNNNNNNNNNNNNNNNNNNNNNgcaatgagAT'


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
