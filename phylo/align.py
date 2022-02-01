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

import collections
import subprocess
import statistics
import sys

from .alignment import Alignment
from .log import log, section_header, explanation
from .misc import count_seqs_in_fasta


def align(args):
    welcome_message()
    samples = find_samples(args.in_dir)
    build_indices(args.in_dir, samples)
    align_all_samples(args.in_dir, args.out_file, samples, args.threads)
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


def find_samples(in_dir):
    samples = sorted(f.name[:-13] for f in in_dir.glob('*_pieces.fasta') if f.is_file())
    if not samples:
        sys.exit(f'\nError: no *_pieces.fasta files found in {in_dir}')
    for s in samples:
        pieces_fasta = in_dir / (s + '_pieces.fasta')
        sequence_fasta = in_dir / (s + '_full.fasta')
        if not pieces_fasta.is_file():
            sys.exit(f'\nError: {pieces_fasta} file not found')
        if not sequence_fasta.is_file():
            sys.exit(f'\nError: {sequence_fasta} file not found')
    log(f'Found {len(samples):,} samples in {in_dir.resolve()}')
    log()
    return samples


def build_indices(in_dir, samples):
    section_header('Building alignment indices')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')
    # TODO: do this in parallel in a thread pool?
    log(f'0 / {len(samples)}', end='')
    for i, s in enumerate(samples):
        if not index_exists(in_dir, s):
            sequence_fasta = (in_dir / (s + '_full.fasta')).resolve()
            index_file = (in_dir / (s + '_full.mmi')).resolve()
            # TODO: explore different indexing options (e.g. -k and -w) to see how they affect
            #       the results.
            command = ['minimap2', '-k15', '-w10',
                       '-d', index_file, sequence_fasta]
            p = subprocess.run(command, capture_output=True, text=True)
            if p.returncode != 0:
                sys.exit(f'\nError: minimap2 failed to index sample {s}:\n{p.stderr}')
        log(f'\r{i+1} / {len(samples)}', end='')
    log('\n')


def index_exists(in_dir, sample):
    index = in_dir / (sample + '_full.mmi')
    return index.is_file() and index.stat().st_size > 0


def align_all_samples(in_dir, out_filename, samples, threads):
    section_header('Aligning pairwise combinations')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')
    with open(out_filename, 'wt') as out_file:
        for a in samples:
            for b in samples:
                if a != b:
                    align_sample_pair(in_dir, out_file, a, b, threads)


def align_sample_pair(in_dir, out_file, a, b, threads):
    pieces_fasta = in_dir / (a + '_pieces.fasta')
    sequence_index = in_dir / (b + '_full.mmi')
    piece_count = count_seqs_in_fasta(pieces_fasta)
    log(f'Aligning {a} pieces to {b} assembly:')
    # TODO: explore different alignment options (e.g. the things set by -x sr) to see how they
    #       affect the results.
    command = ['minimap2', '-a', '-t', str(threads), '--eqx', '-x', 'sr',
               str(sequence_index.resolve()), str(pieces_fasta.resolve())]
    p = subprocess.run(command, capture_output=True, text=True)

    alignments = [Alignment(line) for line in p.stdout.splitlines() if not line.startswith('@')]
    alignments = [a for a in alignments if a.is_fully_aligned()]
    alignments = get_best_alignment_per_read(alignments)
    aligned_fraction = len(alignments) / piece_count
    log(f'  {len(alignments):,} / {piece_count:,} pieces fully aligned '
        f'({100 * aligned_fraction:.1f}%)')

    match_counts = [a.get_match_count() for a in alignments]
    piece_len = int(statistics.median(a.read_length for a in alignments))
    median_identity = 100.0 * statistics.median(match_counts) / piece_len
    mean_identity = 100.0 * statistics.mean(match_counts) / piece_len

    log(f'  median identity: {median_identity:6.2f}%')
    log(f'  mean identity:   {mean_identity:6.2f}%')
    log()

    distances = [piece_len - m for m in match_counts]
    distance_counts = collections.Counter(distances)

    out_file.write(f'{a}\t{b}\t{piece_len}\t{aligned_fraction:.8f}')
    for i in range(max(distances) + 1):
        if distance_counts[i] == 0:
            out_file.write('\t0')
        else:
            probability_mass = distance_counts[i] / piece_count
            out_file.write(f'\t{probability_mass:.8f}')
    out_file.write('\n')


def get_best_alignment_per_read(alignments):
    """
    Returns a list of alignments where there is only one alignment per read. 'Best' is defined using
    the match count.
    """
    alignments_by_read = {}
    for a in alignments:
        if a.read_name not in alignments_by_read:
            alignments_by_read[a.read_name] = a
        elif a.get_match_count() > alignments_by_read[a.read_name].get_match_count():
            alignments_by_read[a.read_name] = a
    return list(alignments_by_read.values())

