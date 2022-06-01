"""
This module contains some miscellaneous functions related to reading and parsing the TSV file made
by Verticall pairwise.

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

import sys


def check_first_two_columns(header_parts, filename):
    """
    This function checks that the first two columns of the header are 'assembly_a' and 'assembly_b'.
    If not, the given file is probably not a tsv from Verticall pairwise.
    """
    if header_parts[0] != 'assembly_a':
        sys.exit(f'Error: first column in {filename} is not labelled "assembly_a" - is the file '
                 f'formatted correctly?')
    if header_parts[1] != 'assembly_b':
        sys.exit(f'Error: second column in {filename} is not labelled "assembly_b" - is the file '
                 f'formatted correctly?')


def check_specific_column(column, header_parts, filename):
    if column not in header_parts:
        sys.exit(f'Error: no column named "{column}" found in {filename} - is the file formatted '
                 f'correctly?')


def get_column_index(header_parts, column, filename):
    check_first_two_columns(header_parts, filename)
    check_specific_column(column, header_parts, filename)
    for i, header in enumerate(header_parts):
        if column == header:
            return i


def check_header_for_assembly_a_regions(header_parts, filename):
    check_first_two_columns(header_parts, filename)
    check_specific_column('assembly_a_vertical_regions', header_parts, filename)
    check_specific_column('assembly_a_horizontal_regions', header_parts, filename)
    check_specific_column('assembly_a_unaligned_regions', header_parts, filename)


def split_region_str(region):
    """
    Input:  a contig name/region string, e.g. "contig_1:123-456"
    Output: a tuple of the contig name and start/end positions, e.g. (contig_1, 123, 456)
    """
    try:
        name, contig_range = region.split(':')
        start, end = contig_range.split('-')
        return name, int(start), int(end)
    except ValueError:
        sys.exit(f'Error: data is not correctly formatted: {region}')


def get_start_end(region):
    """
    Does the same thing as split_region_str, but doesn't include the contig name in the result.
    Input:  a contig name/region string, e.g. "contig_1:123-456"
    Output: a tuple of ints with the start/end positions, e.g. (123, 456)
    """
    try:
        contig_name, range_str = region.split(':')
        start_pos, end_pos = range_str.split('-')
        return int(start_pos), int(end_pos)
    except ValueError:
        sys.exit(f'Error: data is not correctly formatted: {region}')
