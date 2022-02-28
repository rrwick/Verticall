"""
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

import collections
import re
import subprocess
import sys

from .intrange import IntRange
from .log import log, section_header, explanation
from .misc import get_fasta_size, get_n50


def build_indices(in_dir, assemblies):
    section_header('Building alignment indices')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')
    # TODO: do this in parallel in a process pool?
    log(f'0 / {len(assemblies)}', end='')
    for i, a in enumerate(assemblies):
        sample_name, assembly_filename = a
        if not index_exists(in_dir, sample_name):
            index_file = (in_dir / (sample_name + '.mmi')).resolve()
            # TODO: explore different indexing options (e.g. -k and -w) to see how they affect
            #       the results.
            command = ['minimap2', '-k15', '-w10', '-d', index_file, assembly_filename]
            p = subprocess.run(command, capture_output=True, text=True)
            if p.returncode != 0:
                sys.exit(f'\nError: minimap2 failed to index sample {sample_name}:\n{p.stderr}')
        log(f'\r{i+1} / {len(assemblies)}', end='')
    log('\n')


def index_exists(in_dir, sample_name):
    index = in_dir / (sample_name + '.mmi')
    # TODO: is there a way to check for malformed indices?
    return index.is_file() and index.stat().st_size > 0


def align_sample_pair(args, assembly_filename_a, sample_name_b):
    sequence_index = args.in_dir / (sample_name_b + '.mmi')
    # TODO: explore different alignment options (e.g. the things set by -x asm20) to see how they
    #       affect the results.
    # TODO: make minimap2 alignment options settable via an option
    command = ['minimap2', '-c', '-t', '1', '--eqx', '-x', 'asm20',
               str(sequence_index.resolve()), str(assembly_filename_a.resolve())]
    p = subprocess.run(command, capture_output=True, text=True)

    alignments = [Alignment(line) for line in p.stdout.splitlines() if not line.startswith('@')]
    alignments = cull_redundant_alignments(alignments, args.allowed_overlap)

    if not alignments:
        return [], ['  no alignments found']

    return alignments, [f'  {len(alignments)} alignments']


def get_distribution(args, alignments, assembly_filename_a):
    """
    Uses the alignments to build a distance distribution.
    """
    if args.ignore_indels:
        all_cigars = [remove_indels(a.expanded_cigar) for a in alignments]
    else:
        all_cigars = [compress_indels(a.expanded_cigar) for a in alignments]

    n50_alignment_length = get_n50(len(c) for c in all_cigars)
    log_text = [f'  N50 alignment length: {n50_alignment_length}']

    query_coverage = get_query_coverage(alignments, assembly_filename_a)
    mean_identity = 1.0 - (sum(get_difference_count(c) for c in all_cigars) /
                           sum(len(c) for c in all_cigars))

    log_text.append(f'  aligned fraction: {100.0 * query_coverage:.2f}%')
    log_text.append(f'  mean identity: {100.0 * mean_identity:.2f}%')

    window_size, window_step = choose_window_size_and_step(all_cigars, args.window_count)
    all_cigars = [c for c in all_cigars if len(c) >= window_size]

    distances, max_difference_count = get_distances(all_cigars, window_size, window_step)
    distance_counts = collections.Counter(distances)
    log_text.append(f'  distances sampled from {len(distances)} x {window_size} bp windows')

    return distance_counts, query_coverage, log_text


def cull_redundant_alignments(alignments, allowed_overlap):
    alignments = sorted(alignments, key=lambda x: x.alignment_score, reverse=True)
    alignments_by_contig = collections.defaultdict(list)
    alignments_no_redundancy = []
    for a in alignments:
        if not any(a.overlaps(b, allowed_overlap) for b in alignments_by_contig[a.query_name]):
            alignments_no_redundancy.append(a)
        alignments_by_contig[a.query_name].append(a)
    return alignments_no_redundancy


def get_query_coverage(alignments, assembly_filename):
    assembly_size = get_fasta_size(assembly_filename)
    ranges_by_contig = {}
    for a in alignments:
        if a.query_name not in ranges_by_contig:
            ranges_by_contig[a.query_name] = IntRange()
        ranges_by_contig[a.query_name].add_range(a.query_start, a.query_end)
    aligned_bases = sum(r.total_length() for r in ranges_by_contig.values())
    assert aligned_bases <= assembly_size
    return aligned_bases / assembly_size


def get_distances(all_cigars, window_size, window_step):
    distances, max_difference_count = [], 0
    for cigar in all_cigars:
        # TODO: trim CIGAR so windows fit in the middle?
        start, end = 0, window_size
        while end <= len(cigar):
            cigar_window = cigar[start:end]
            assert len(cigar_window) == window_size
            difference_count = get_difference_count(cigar_window)
            distances.append(difference_count / window_size)
            max_difference_count = max(max_difference_count, difference_count)
            start += window_step
            end += window_step
    return distances, max_difference_count


def remove_indels(cigar):
    """
    Removes all indels from a CIGAR. For example:
    in:  ===X=IIII==XX==DDDD==
    out: ===X===XX====
    """
    return re.sub(r'I+', '', re.sub(r'D+', '', cigar))


def compress_indels(cigar):
    """
    Compresses runs of indels into size-1 indels. For example:
    in:  ===X=IIII==XX==DDDD==
    out: ===X=I==XX==D==
    """
    return re.sub(r'I+', 'I', re.sub(r'D+', 'D', cigar))


def get_difference_count(cigar):
    """
    Returns the number of mismatches and indels in the CIGAR.
    """
    return cigar.count('X') + cigar.count('I') + cigar.count('D')


def choose_window_size_and_step(cigars, target_window_count):
    """
    This function chooses an appropriate window size and step for the given CIGARs. It tries to
    balance larger windows, which give higher-resolution identity samples, especially with
    closely-related assemblies, and smaller windows, which allow for more identity samples.
    """
    window_step = 1000
    while window_step > 1:
        window_size = window_step * 100
        if get_window_count(cigars, window_size, window_step) > target_window_count:
            return window_size, window_step
        window_step -= 1
    return window_step * 100, window_step


def get_window_count(cigars, window_size, window_step):
    """
    For a given window size, window step and set of CIGARs, this function returns how many windows
    there will be in total.
    """
    count = 0
    for cigar in cigars:
        cigar_len = len(cigar)
        if cigar_len < window_size:
            continue
        cigar_len -= window_size
        count += 1
        count += cigar_len // window_step
    return count


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
