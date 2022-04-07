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

import verticall.repair


def test_split_seq_on_ambiguous():
    assert verticall.repair.split_seq_on_ambiguous('ACGACTACGACTACG') == ['ACGACTACGACTACG']
    assert verticall.repair.split_seq_on_ambiguous('ACGACTNCGACTACG') == ['ACGACT', 'CGACTACG']
    assert verticall.repair.split_seq_on_ambiguous('ACGACNNNNNCTACG') == ['ACGAC', 'CTACG']
    assert verticall.repair.split_seq_on_ambiguous('ACGACNACGNCTNNG') == ['ACGAC', 'ACG', 'CT', 'G']
    assert verticall.repair.split_seq_on_ambiguous('NNNNNNNNNNNNNNN') == []
    assert verticall.repair.split_seq_on_ambiguous('ACG?CTAQGXXTACG') == ['ACG', 'CTA', 'G', 'TACG']
    assert verticall.repair.split_seq_on_ambiguous('RRRACTACGANNNNN') == ['ACTACGA']


def test_make_names_unique():
    assert verticall.repair.make_names_unique(['A', 'B', 'C']) == ['A', 'B', 'C']
    assert verticall.repair.make_names_unique(['A', 'B', 'B']) == ['A', 'B_1', 'B_2']
    assert verticall.repair.make_names_unique(['B', 'B', 'B']) == ['B_1', 'B_2', 'B_3']
    assert verticall.repair.make_names_unique(['A', 'B', 'A']) == ['A_1', 'B', 'A_2']
    assert verticall.repair.make_names_unique(['A', 'A', 'A_1']) == ['A_1_1', 'A_2', 'A_1_2']
