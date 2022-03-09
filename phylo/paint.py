"""
Copyright 2022 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/XXXXXXXXX

This module contains code related to 'painting' an assembly using the alignments and the high/low
identity thresholds.

This file is part of XXXXXXXXX. XXXXXXXXX is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. XXXXXXXXX is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with XXXXXXXXX.
If not, see <https://www.gnu.org/licenses/>.
"""

import enum

from .distance import get_vertical_horizontal_distributions
from .intrange import IntRange
from .misc import iterate_fasta


class AlignmentRole(enum.Enum):
    QUERY = 0
    TARGET = 1


class Paint(enum.Enum):
    UNALIGNED = 0
    VERTICAL = 1
    HORIZONTAL = 2
    AMBIGUOUS = 3

    def __repr__(self):
        if self == Paint.UNALIGNED:
            return 'U'
        elif self == Paint.VERTICAL:
            return 'V'
        elif self == Paint.HORIZONTAL:
            return 'H'
        elif self == Paint.AMBIGUOUS:
            return '?'
        else:
            assert False


def paint_alignments(alignments, thresholds):
    log_text = ['  painting alignments:']
    for a in alignments:
        a.paint_sliding_windows(thresholds)
    vertical_masses, horizontal_masses = get_vertical_horizontal_distributions(alignments)
    total_vertical_mass = sum(vertical_masses)
    total_horizontal_mass = sum(horizontal_masses)
    log_text.append(f'    vertical inheritance: {100.0 * total_vertical_mass:.2f}%')
    log_text.append(f'    horizontal inheritance: {100.0 * total_horizontal_mass:.2f}%')
    return vertical_masses, horizontal_masses, log_text


def paint_assemblies(name_a, name_b, filename_a, filename_b, alignments):
    painted_a = PaintedAssembly(filename_a)
    painted_b = PaintedAssembly(filename_b)

    for a in alignments:
        painted_a.add_alignment(a, AlignmentRole.QUERY)
        painted_b.add_alignment(a, AlignmentRole.TARGET)

    log_text = get_painting_log_text(name_a, painted_a)
    log_text += get_painting_log_text(name_b, painted_b)

    return painted_a, painted_b, log_text


def get_painting_log_text(name, painted):
    vertical, horizontal, unaligned = painted.get_fractions()
    log_text = [f'  painting {name}:',
                f'    vertical:   {100.0 * vertical:.2f}%',
                f'    horizontal: {100.0 * horizontal:.2f}%',
                f'    unaligned:  {100.0 * unaligned:.2f}%']
    return log_text


class PaintedAssembly(object):

    def __init__(self, fasta_filename):
        self.contigs = {}
        for name, seq in iterate_fasta(fasta_filename):
            self.contigs[name] = PaintedContig(seq)

    def add_alignment(self, a, role):
        name = a.query_name if role == AlignmentRole.QUERY else a.target_name
        self.contigs[name].add_alignment(a, role)

    def get_max_differences(self):
        """
        Returns the largest number of differences in all of the alignment's sliding windows.
        """
        if len(self.contigs) == 0:
            return 0
        else:
            return max(c.get_max_differences() for c in self.contigs.values())

    def get_fractions(self):
        """
        Returns the vertical, horizontal and unaligned fractions of the assembly.
        """
        total, vertical, horizontal = 0, 0, 0
        for c in self.contigs.values():
            total += c.length
            for start, end in c.get_vertical_blocks():
                vertical += (end - start)
            for start, end in c.get_horizontal_blocks():
                horizontal += (end - start)
        unaligned = total - vertical - horizontal
        return vertical/total, horizontal/total, unaligned/total


class PaintedContig(object):

    def __init__(self, seq):
        self.length = len(seq)
        self.paint = [Paint.UNALIGNED] * self.length
        self.alignment_points = []

    def add_alignment(self, a, role):
        cigar_to_seq = a.cigar_to_query if role == AlignmentRole.QUERY else a.cigar_to_target
        points = []
        vertical_ranges = IntRange()
        horizontal_ranges = IntRange()
        for i, window in enumerate(a.windows_no_overlap):
            differences = a.window_differences[i]
            classification = a.window_classifications[i]
            a_start, a_end = window
            seq_start, seq_end = cigar_to_seq[a_start], cigar_to_seq[a_end-1]+1
            seq_centre = (seq_start + seq_end) / 2
            points.append((seq_centre, differences))

            if classification == Paint.VERTICAL:
                vertical_ranges.add_range(seq_start, seq_end)
            elif classification == Paint.HORIZONTAL:
                horizontal_ranges.add_range(seq_start, seq_end)
            else:
                assert False

        for start, end in horizontal_ranges.ranges:
            for i in range(start, end):
                self.paint_position(i, Paint.HORIZONTAL)
        for start, end in vertical_ranges.ranges:
            for i in range(start, end):
                self.paint_position(i, Paint.VERTICAL)

        self.alignment_points.append(points)

    def get_max_differences(self):
        max_differences = 0
        for points in self.alignment_points:
            differences = [p[1] for p in points]
            if differences:
                max_differences = max(max_differences, max(differences))
        return max_differences

    def paint_position(self, i, classification):
        """
        Paints a single position of the contig. Both VERTICAL and HORIZONTAL paint over UNALIGNED,
        and VERTICAL paints over HORIZONTAL. I.e. VERTICAL takes precedence, then HORIZONTAL, then
        UNALIGNED.
        """
        assert classification == Paint.VERTICAL or classification == Paint.HORIZONTAL
        if self.paint[i] != Paint.VERTICAL:
            self.paint[i] = classification

    def get_vertical_blocks(self):
        """
        Returns a list of all ranges of the contig which have been painted as vertical.
        """
        return get_blocks(self.paint, Paint.VERTICAL)

    def get_horizontal_blocks(self):
        """
        Returns a list of all ranges of the contig which have been painted as horizontal.
        """
        return get_blocks(self.paint, Paint.HORIZONTAL)


def get_blocks(paint, classification):
    blocks = []
    start, end = None, None
    for i, c in enumerate(paint):
        if c == classification:
            if start is None:
                start = i
            end = i + 1
        elif start is not None:
            blocks.append((start, end))
            start, end = None, None
    if start is not None:
        blocks.append((start, end))
    return blocks
