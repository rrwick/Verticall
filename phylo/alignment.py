"""
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

import re
import sys


class Alignment(object):

    def __init__(self, paf_line):
        parts = paf_line.strip().split('\t')
        if len(parts) < 11:
            sys.exit('\nError: alignment file does not seem to be in SAM format')

        self.query_name = parts[0]
        self.query_length = int(parts[1])
        self.query_start = int(parts[2])
        self.query_end = int(parts[3])
        self.strand = parts[4]

        self.target_name = parts[5]
        self.target_length = int(parts[6])
        self.target_start = int(parts[7])
        self.target_end = int(parts[8])

        self.matches = int(parts[9])
        self.alignment_length = int(parts[10])
        self.percent_identity = 100.0 * self.matches / self.alignment_length

        self.cigar, self.alignment_score = None, None
        for part in parts:
            if part.startswith('cg:Z:'):
                self.cigar = part[5:]
            if part.startswith('AS:i:'):
                self.alignment_score = int(part[5:])
        self.expanded_cigar = get_expanded_cigar(self.cigar)

    def query_covered_bases(self):
        return self.query_end - self.query_start

    def __repr__(self):
        return self.query_name + ':' + str(self.query_start) + '-' + str(self.query_end) + \
               '(' + self.strand + '), ' + \
               self.target_name + ':' + str(self.target_start) + '-' + str(self.target_end) + \
               ' (' + ('%.3f' % self.percent_identity) + '%)'

    def overlaps(self, other, allowed_overlap):
        """
        Tests whether this alignment overlaps with the other alignment in the query sequence. A bit
        of overlap can be allowed using the allowed_overlap parameter.
        """
        if self.query_name != other.query_name:
            return False
        this_start = self.query_start + allowed_overlap
        this_end = self.query_end - allowed_overlap
        if this_start >= this_end:
            return False
        return (this_start < other.query_end) and (other.query_start < this_end)


def get_expanded_cigar(cigar):
    expanded_cigar = []
    cigar_parts = re.findall(r'\d+[IDX=]', cigar)
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        expanded_cigar.append(letter * size)
    return ''.join(expanded_cigar)
