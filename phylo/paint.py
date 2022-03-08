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

from .alignment import remove_indels, compress_indels, swap_insertions_and_deletions, \
    cigar_to_contig_pos, get_difference_count
from .distance import get_vertical_horizontal_distributions
from .misc import iterate_fasta


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


def paint_assemblies(args, name_a, name_b, filename_a, filename_b, alignments, window_size,
                     thresholds):
    log_text = [f'  painting {name_a}']
    painted_a = PaintedAssembly(filename_a)
    painted_b = PaintedAssembly(filename_b)

    for a in alignments:
        painted_a.add_alignment(a, 'query', window_size, args.ignore_indels)
        painted_b.add_alignment(a, 'target', window_size, args.ignore_indels)

    # TODO: finalise the painting
    #   * assign a value to each position of each contig:
    #     * alignment identity (averaged over windows) if aligned
    #     * 'None' if not aligned
    #   * make a simplified collection of points for plotting (per aligned region)
    #   * define recombinant vs non-recombinant regions
    #     * intermediate positions surrounded by non-recombinant positions are non-recombinant
    #     * intermediate positions surrounded by recombinant positions are recombinant
    #     * intermediate positions surrounded by both (recombinant on one side, non-recombinant on
    #       the other) are recombinant (i.e. erring on the side of calling recombination).

    return painted_a, painted_b, log_text


class PaintedAssembly(object):

    def __init__(self, fasta_filename):
        self.contigs = {}
        for name, seq in iterate_fasta(fasta_filename):
            self.contigs[name] = PaintedContig(seq)

    def add_alignment(self, a, assembly_status, window_size, ignore_indels):
        if assembly_status == 'query':
            name, start, end = a.query_name, a.query_start, a.query_end
            cigar = a.expanded_cigar
        elif assembly_status == 'target':
            name, start, end = a.target_name, a.target_start, a.target_end
            cigar = swap_insertions_and_deletions(a.expanded_cigar)
        else:
            assert False
        self.contigs[name].add_alignment(start, end, cigar, window_size, ignore_indels)

    def finalise(self):
        for c in self.contigs.values():
            c.finalise()

    def get_max_differences(self):
        if len(self.contigs) == 0:
            return 0
        else:
            return max(c.get_max_differences() for c in self.contigs.values())


class PaintedContig(object):

    def __init__(self, seq):
        self.length = len(seq)
        self.differences = [0] * self.length  # difference count per aligned contig position
        self.window_differences = []   # (contig start, contig end, difference count)

    def add_alignment(self, a_start, a_end, cigar, window_size, ignore_indels):
        assert window_size % 100 == 0
        window_step = window_size // 100
        cigar_to_contig = cigar_to_contig_pos(cigar, a_start, a_end)
        if ignore_indels:
            cigar, cigar_to_contig = remove_indels(cigar, cigar_to_contig)
        else:
            cigar, cigar_to_contig = compress_indels(cigar, cigar_to_contig)

        # TODO: rework this loop so we sample to the very end of the alignment?
        start, end = 0, window_size
        while end <= len(cigar):
            cigar_window = cigar[start:end]
            assert len(cigar_window) == window_size
            difference_count = get_difference_count(cigar_window)
            self.window_differences.append((cigar_to_contig[start], cigar_to_contig[end-1],
                                            difference_count))
            start += window_step
            end += window_step

    def finalise(self):
        counts = [None] * self.length
        for contig_range, differences in self.window_differences:
            start, end = contig_range
            # TODO
            # TODO
            # TODO
            # TODO
            # TODO

    def get_max_differences(self):
        if len(self.window_differences) == 0:
            return 0
        else:
            return max(d[2] for d in self.window_differences)
