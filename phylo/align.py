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
from multiprocessing.pool import ThreadPool
import re
import subprocess
import sys

from .alignment import Alignment
from .intrange import IntRange
from .log import log, section_header, explanation
from .misc import get_fasta_size


def align(args):
    welcome_message()
    assemblies = find_assemblies(args.in_dir)
    build_indices(args.in_dir, assemblies)
    align_all_samples(args.in_dir, args.out_file, assemblies, args.threads, args.allowed_overlap,
                      args.window_size, args.window_step, args.ignore_indels)
    finished_message()


def welcome_message():
    section_header('Starting XXXXXXXXX align')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')


def finished_message():
    section_header('Finished!')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')


def find_assemblies(in_dir):
    """
    Returns assemblies in a (sample_name, filename) tuple.
    """
    def find_assemblies_with_extension(extension, all_assemblies):
        extension_len = len(extension) + 1  # plus one for the dot
        fasta_assemblies = sorted(f for f in in_dir.glob('*.' + extension) if f.is_file())
        for a in fasta_assemblies:
            sample_name = a.name[:-extension_len]
            if sample_name in all_assemblies:
                sys.exit(f'\nError: duplicate sample name {sample_name}')
            all_assemblies[sample_name] = a

    assemblies = {}
    find_assemblies_with_extension('fasta', assemblies)
    find_assemblies_with_extension('fasta.gz', assemblies)
    find_assemblies_with_extension('fna', assemblies)
    find_assemblies_with_extension('fna.gz', assemblies)
    assemblies = sorted(assemblies.items())

    log(f'Found {len(assemblies):,} samples in {in_dir.resolve()}')
    log()
    return assemblies


def build_indices(in_dir, assemblies):
    section_header('Building alignment indices')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')
    # TODO: do this in parallel in a thread pool?
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
    return index.is_file() and index.stat().st_size > 0


def align_all_samples(in_dir, out_filename, assemblies, threads, allowed_overlap, window_size,
                      window_step, ignore_indels):
    section_header('Aligning pairwise combinations')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')

    arg_list = []
    for sample_name_a, assembly_filename_a in assemblies:
        for sample_name_b, assembly_filename_b in assemblies:
            if sample_name_a != sample_name_b:
                arg_list.append((in_dir, sample_name_a, sample_name_b, assembly_filename_a,
                                 allowed_overlap, window_size, window_step, ignore_indels))

    with open(out_filename, 'wt') as out_file:
        out_file.write('#sample_a\tsample_b\twindow_size\talignment_coverage\tprobability_masses\n')

        # If only using a single thread, do the alignment in a simple loop (easier for debugging).
        if threads == 1:
            for a in arg_list:
                output_line, log_text = align_sample_pair(a)
                out_file.write(output_line)
                log('\n'.join(log_text), end='\n\n')

        # If only using multiple threads, do the alignments in a thread pool.
        else:
            with ThreadPool(processes=threads) as pool:
                for output_line, log_text in pool.imap(align_sample_pair, arg_list):
                    out_file.write(output_line)
                    log('\n'.join(log_text), end='\n\n')


def align_sample_pair(all_args):
    """
    Arguments are passes as a single tuple to make this function easier to call via pool.imap.
    """
    in_dir, sample_name_a, sample_name_b, assembly_filename_a, allowed_overlap, window_size, \
        window_step, ignore_indels = all_args

    sequence_index = in_dir / (sample_name_b + '.mmi')
    log_text = [f'Aligning {sample_name_a} to {sample_name_b}:']
    # TODO: explore different alignment options (e.g. the things set by -x asm20) to see how they
    #       affect the results.
    command = ['minimap2', '-c', '-t', '1', '--eqx', '-x', 'asm20',
               str(sequence_index.resolve()), str(assembly_filename_a.resolve())]
    p = subprocess.run(command, capture_output=True, text=True)

    alignments = [Alignment(line) for line in p.stdout.splitlines() if not line.startswith('@')]
    full_alignment_count = len(alignments)
    alignments = [a for a in alignments if a.alignment_length >= window_size]
    alignments = cull_redundant_alignments(alignments, allowed_overlap)
    if ignore_indels:
        all_cigars = [remove_indels(a.expanded_cigar) for a in alignments]
    else:
        all_cigars = [compress_indels(a.expanded_cigar) for a in alignments]
    all_cigars = [c for c in all_cigars if len(c) >= window_size]
    log_text.append(f'  {full_alignment_count} alignments ({len(all_cigars)} after filtering)')

    distances, max_difference_count = get_distances(all_cigars, window_size, window_step)
    distance_counts = collections.Counter(distances)
    log_text.append(f'  distance sampled from {len(distances)} alignment windows')

    query_coverage = get_query_coverage(alignments, assembly_filename_a)
    mean_identity = 1.0 - (sum(get_difference_count(c) for c in all_cigars) /
                           sum(len(c) for c in all_cigars))

    log_text.append(f'  aligned fraction: {100.0 * query_coverage:6.2f}%')
    log_text.append(f'  mean identity:    {100.0 * mean_identity:6.2f}%')

    output_line = [sample_name_a, sample_name_b, str(window_size), f'{query_coverage:.8f}']
    for i in range(max_difference_count + 1):
        d = i / window_size
        if distance_counts[d] == 0:
            output_line.append('0')
        else:
            probability_mass = distance_counts[d] / len(distances)
            output_line.append(f'{probability_mass:.8f}')
    return '\t'.join(output_line) + '\n', log_text


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
