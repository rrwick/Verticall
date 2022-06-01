"""
This module contains code related to 'painting' an assembly using the alignments and the high/low
identity thresholds.

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

import enum

from .distance import get_vertical_horizontal_distributions, get_distance
from .misc import iterate_fasta, get_difference_count


class AlignmentRole(enum.Enum):
    QUERY = 0
    TARGET = 1


def paint_alignments(alignments, thresholds, window_size):
    for a in alignments:
        a.paint_sliding_windows(thresholds)
    vertical_masses, horizontal_masses = get_vertical_horizontal_distributions(alignments)
    total_vertical_mass = sum(vertical_masses)
    total_horizontal_mass = sum(horizontal_masses)
    mean_vert_window_dist = get_distance(vertical_masses, window_size, 'mean')
    median_vert_window_dist = get_distance(vertical_masses, window_size, 'median')
    mean_vert_dist = get_mean_vertical_distance(alignments)
    r_over_m = get_r_over_m(alignments)

    log_text = ['  painting alignments:',
                f'    vertical inheritance:          {100.0 * total_vertical_mass:6.2f}%',
                f'    horizontal inheritance:        {100.0 * total_horizontal_mass:6.2f}%',
                f'    mean vertical window distance:   {mean_vert_window_dist:.9f}',
                f'    median vertical window distance: {median_vert_window_dist:.9f}',
                f'    mean vertical distance:          {mean_vert_dist:.9f}']
    return vertical_masses, horizontal_masses, mean_vert_window_dist, median_vert_window_dist, \
        mean_vert_dist, r_over_m, log_text


def get_mean_vertical_distance(alignments):
    """
    This function uses just the vertically-painted regions to get a mean distance.
    """
    total_size, differences = 0, 0
    for a in alignments:
        for start, end in a.get_vertical_blocks():
            total_size += (end - start)
            differences += get_difference_count(a.simplified_cigar[start:end])
    if total_size == 0:
        return 0.0
    else:
        return differences / total_size


def get_r_over_m(alignments):
    v_differences, h_differences = 0, 0
    for a in alignments:
        for start, end in a.get_vertical_blocks():
            v_differences += get_difference_count(a.simplified_cigar[start:end])
        for start, end in a.get_horizontal_blocks():
            h_differences += get_difference_count(a.simplified_cigar[start:end])
    if v_differences == 0:
        return 0.0
    else:
        return h_differences / v_differences


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
    log_text = [f'  painting {name} contigs:',
                f'    vertical:   {100.0 * vertical:6.2f}%',
                f'    horizontal: {100.0 * horizontal:6.2f}%',
                f'    unaligned:  {100.0 * unaligned:6.2f}%']
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
        if total == 0:
            return 0.0, 0.0, 0.0
        return vertical/total, horizontal/total, unaligned/total

    def get_regions(self):
        """
        Returns a string encoding the vertical, horizontal and unaligned regions of the assembly.
        """
        vertical, horizontal, unaligned = [], [], []
        for name, c in self.contigs.items():
            for start, end in c.get_vertical_blocks():
                vertical.append(f'{name}:{start}-{end}')
            for start, end in c.get_horizontal_blocks():
                horizontal.append(f'{name}:{start}-{end}')
            for start, end in c.get_unaligned_blocks():
                unaligned.append(f'{name}:{start}-{end}')
        return ','.join(vertical), ','.join(horizontal), ','.join(unaligned)


class PaintedContig(object):

    def __init__(self, seq):
        self.length = len(seq)
        self.paint = [0] * self.length  # 0 means unaligned
        self.alignment_points = []
        self.vertical_blocks = None
        self.horizontal_blocks = None
        self.unaligned_blocks = None

    def add_alignment(self, a, role):
        cigar_to_seq = a.cigar_to_query if role == AlignmentRole.QUERY else a.cigar_to_target
        points = []
        vertical_ranges, horizontal_ranges = [], []
        for i, window in enumerate(a.windows_no_overlap):
            differences = a.window_differences[i]
            classification = a.window_classifications[i]
            a_start, a_end = window

            seq_start, seq_end = cigar_to_seq[a_start], cigar_to_seq[a_end-1]
            if seq_end < seq_start:
                seq_start, seq_end = seq_end, seq_start
            seq_end += 1

            seq_centre = (seq_start + seq_end) / 2
            points.append((seq_centre, differences))

            if classification == 1:  # 1 means vertical
                vertical_ranges.append((seq_start, seq_end))
            elif classification == 2:  # 2 means vertical
                horizontal_ranges.append((seq_start, seq_end))
            else:
                assert False

        # Both vertical (1) and horizontal (2) paint over unaligned (3), and vertical paints over
        # horizontal. I.e. vertical takes precedence, then horizontal, then unaligned.
        for start, end in horizontal_ranges:
            for i in range(start, end):
                if self.paint[i] != 1:  # 1 means vertical
                    self.paint[i] = 2  # 2 means horizontal
        for start, end in vertical_ranges:
            for i in range(start, end):
                self.paint[i] = 1  # 1 means vertical

        self.alignment_points.append(points)

    def get_max_differences(self):
        max_differences = 0
        for points in self.alignment_points:
            differences = [p[1] for p in points]
            if differences:
                max_differences = max(max_differences, max(differences))
        return max_differences

    def get_vertical_blocks(self):
        """
        Returns a list of all ranges of the contig which have been painted as vertical.
        """
        if self.vertical_blocks is None:
            self.vertical_blocks = get_blocks(self.paint, 1)  # 1 means vertical
        return self.vertical_blocks

    def get_horizontal_blocks(self):
        """
        Returns a list of all ranges of the contig which have been painted as horizontal.
        """
        if self.horizontal_blocks is None:
            self.horizontal_blocks = get_blocks(self.paint, 2)  # 2 means horizontal
        return self.horizontal_blocks

    def get_unaligned_blocks(self):
        """
        Returns a list of all ranges of the contig which have been painted as unaligned.
        """
        if self.unaligned_blocks is None:
            self.unaligned_blocks = get_blocks(self.paint, 0)  # 0 means unaligned
        return self.unaligned_blocks


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
