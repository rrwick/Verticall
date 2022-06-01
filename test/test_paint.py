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

import verticall.paint


def test_get_blocks():
    paint = [0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0]
    assert verticall.paint.get_blocks(paint, 0) == [(0, 1), (4, 7), (8, 9), (11, 12)]
    assert verticall.paint.get_blocks(paint, 1) == [(1, 4), (7, 8), (9, 11)]
    assert verticall.paint.get_blocks(paint, 2) == []

    paint = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]
    assert verticall.paint.get_blocks(paint, 0) == [(0, 6)]
    assert verticall.paint.get_blocks(paint, 1) == [(6, 12)]

    paint = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    assert verticall.paint.get_blocks(paint, 0) == [(0, 12)]
    assert verticall.paint.get_blocks(paint, 1) == []

    paint = []
    assert verticall.paint.get_blocks(paint, 0) == []
    assert verticall.paint.get_blocks(paint, 1) == []

    paint = ['A', 'A', 'B', 'C', 'C', 'A']
    assert verticall.paint.get_blocks(paint, 'A') == [(0, 2), (5, 6)]
    assert verticall.paint.get_blocks(paint, 'B') == [(2, 3)]
    assert verticall.paint.get_blocks(paint, 'C') == [(3, 5)]
