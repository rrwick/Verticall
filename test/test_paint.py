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
If not, see <https://www.gnu.org/licenses/>.
"""

import phylo.paint


def test_get_blocks():
    paint = [0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0]
    assert phylo.paint.get_blocks(paint, 0) == [(0, 1), (4, 7), (8, 9), (11, 12)]
    assert phylo.paint.get_blocks(paint, 1) == [(1, 4), (7, 8), (9, 11)]
