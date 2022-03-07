"""
Copyright 2022 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/XXXXXXXXX

This module contains code related to assembly-vs-assembly alignment via minimap2.

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


def build_indices(args, assemblies):
    section_header('Building alignment indices')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')
    index_options = args.index_options.split()
    if not args.verbose:
        log(f'0 / {len(assemblies)}', end='')
    for i, a in enumerate(assemblies):
        sample_name, assembly_filename = a
        if not index_exists(args, sample_name):
            index_file = (args.in_dir / (sample_name + '.mmi')).resolve()
            command = ['minimap2']
            command += index_options
            command += ['-d', index_file, assembly_filename]
            if args.verbose:
                log(' '.join(str(x) for x in command))
            p = subprocess.run(command, capture_output=True, text=True)
            if p.returncode != 0:
                sys.exit(f'\nError: minimap2 failed to index sample {sample_name}:\n{p.stderr}')
        if not args.verbose:
            log(f'\r{i+1} / {len(assemblies)}', end='')
    if not args.verbose:
        log()
    log()


def index_exists(args, sample_name):
    index = args.in_dir / (sample_name + '.mmi')
    # TODO: is there a way to check for malformed indices?
    exists = (index.is_file() and index.stat().st_size > 0)
    if exists and args.verbose:
        log(f'{index} already exists')
    return exists


def align_sample_pair(args, assembly_filename_a, sample_name_b):
    log_text = []
    sequence_index = args.in_dir / (sample_name_b + '.mmi')
    command = ['minimap2', '-c', '-t', '1', '--eqx']
    command += args.align_options.split()
    command += [str(sequence_index.resolve()), str(assembly_filename_a.resolve())]
    if args.verbose:
        log_text.append('  ' + ' '.join(str(x) for x in command))
    p = subprocess.run(command, capture_output=True, text=True)

    ignore_indels = True if args.ignore_indels else False
    alignments = [Alignment(line, ignore_indels) for line in p.stdout.splitlines()
                  if not line.startswith('@')]
    alignments = cull_redundant_alignments(alignments, args.allowed_overlap)

    n50_alignment_length = get_n50(len(a.expanded_cigar) for a in alignments)
    aligned_frac = get_query_coverage(alignments, assembly_filename_a)

    if not alignments:
        log_text.append('  no alignments found')
    else:
        log_text.append(f'  {len(alignments)} alignments:')
        log_text.append(f'    N50 alignment length: {n50_alignment_length}')
        log_text.append(f'    aligned fraction: {100.0 * aligned_frac:.2f}%')
    return alignments, aligned_frac, log_text


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


class Alignment(object):

    def __init__(self, paf_line, ignore_indels=False):
        self.query_name, self.query_length, self.query_start, self.query_end, self.strand, \
            self.target_name, self.target_length, self.target_start, self.target_end, \
            self.matches, self.alignment_length, self.percent_identity, self.cigar, \
            self.alignment_score = self.read_paf_columns(paf_line)
        self.expanded_cigar, self.simplified_cigar, self.cigar_to_query, self.cigar_to_target = \
            self.get_cigars(ignore_indels)

    @staticmethod
    def read_paf_columns(paf_line):
        parts = paf_line.strip().split('\t')
        if len(parts) < 11:
            sys.exit('\nError: alignment file does not seem to be in PAF format')
        query_name = parts[0]
        query_length = int(parts[1])
        query_start = int(parts[2])
        query_end = int(parts[3])
        strand = parts[4]
        target_name = parts[5]
        target_length = int(parts[6])
        target_start = int(parts[7])
        target_end = int(parts[8])
        matches = int(parts[9])
        alignment_length = int(parts[10])
        percent_identity = 100.0 * matches / alignment_length
        cigar, alignment_score = None, None
        for part in parts:
            if part.startswith('cg:Z:'):
                cigar = part[5:]
            if part.startswith('AS:i:'):
                alignment_score = int(part[5:])
        return query_name, query_length, query_start, query_end, strand, target_name,\
            target_length, target_start, target_end, matches, alignment_length, percent_identity,\
            cigar, alignment_score

    def get_cigars(self, ignore_indels):
        expanded_cigar = get_expanded_cigar(self.cigar)
        cigar_to_query = cigar_to_contig_pos(expanded_cigar, self.query_start, self.query_end)
        flipped_cigar = swap_insertions_and_deletions(expanded_cigar)
        cigar_to_target = cigar_to_contig_pos(flipped_cigar, self.target_start, self.target_end)

        # Compress/remove indels from the CIGAR to make a simplified CIGAR over which the sliding
        # window will operate.
        indel_func = remove_indels if ignore_indels else compress_indels
        simplified_cigar, cigar_to_query = indel_func(expanded_cigar, cigar_to_contig=cigar_to_query)
        _, cigar_to_target = indel_func(expanded_cigar, cigar_to_contig=cigar_to_target)
        assert len(simplified_cigar) == len(cigar_to_query) == len(cigar_to_target)

        return expanded_cigar, simplified_cigar, cigar_to_query, cigar_to_target

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


def cigar_to_contig_pos(cigar, start, end):
    """
    Takes an expanded CIGAR and returns a list of the same length, where the list contains
    contig positions corresponding to each CIGAR position.

    Indels are interpreted from the contig's point of view. I.e. an insertion means the contig has
    an extra base and a deletion means the contig is missing a base, so insertions 'consume'
    sequence positions while deletions do not.
    """
    cigar_to_contig = []
    i = start
    for c in cigar:
        cigar_to_contig.append(i)
        if c == '=' or c == 'X' or c == 'I':
            i += 1
    assert i == end
    return cigar_to_contig


def remove_indels(cigar, cigar_to_contig=None):
    """
    Removes all indels from an expanded CIGAR. For example:
    in:  ===X=IIII==XX==DDDD==
    out: ===X===XX====

    Can optionally take a list of cigar-to-contig positions, in which case it will also modify that
    to match the returned CIGAR.
    """
    if cigar_to_contig is None:
        return re.sub(r'I+', '', re.sub(r'D+', '', cigar))
    assert len(cigar) == len(cigar_to_contig)
    new_cigar, new_cigar_to_contig = [], []
    for c, i in zip(cigar, cigar_to_contig):
        if c == '=' or c == 'X':
            new_cigar.append(c)
            new_cigar_to_contig.append(i)
    return ''.join(new_cigar), new_cigar_to_contig


def compress_indels(cigar, cigar_to_contig=None):
    """
    Compresses runs of indels into size-1 indels in an expanded CIGAR. For example:
    in:  ===X=IIII==XX==DDDD==
    out: ===X=I==XX==D==

    Can optionally take a list of cigar-to-contig positions, in which case it will also modify that
    to match the returned CIGAR.
    """
    if cigar_to_contig is None:
        return re.sub(r'I+', 'I', re.sub(r'D+', 'D', cigar))
    assert len(cigar) == len(cigar_to_contig)
    new_cigar, new_cigar_to_contig = [], []
    for c, i in zip(cigar, cigar_to_contig):
        if c == 'D' and len(new_cigar) > 0 and new_cigar[-1] == 'D':
            new_cigar.pop()
            new_cigar_to_contig.pop()
        if c == 'I' and len(new_cigar) > 0 and new_cigar[-1] == 'I':
            new_cigar.pop()
            new_cigar_to_contig.pop()
        new_cigar.append(c)
        new_cigar_to_contig.append(i)
    return ''.join(new_cigar), new_cigar_to_contig


def swap_insertions_and_deletions(cigar):
    """
    Swaps I and D characters in an expanded CIGAR.
    """
    return cigar.replace('I', 'd').replace('D', 'I').replace('d', 'D')
