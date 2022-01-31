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

import pathlib
import subprocess
import sys

from .log import log, section_header, explanation


def align(args):
    welcome_message()
    samples = find_samples(args.in_dir)
    build_indices(args.in_dir, samples)
    align_all_samples(args.in_dir, samples)
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
            sequence_fasta = in_dir / (s + '_full.fasta')
            command = ['bwa', 'index', str(sequence_fasta.resolve())]
            p = subprocess.run(command, capture_output=True, text=True)
            if p.returncode != 0:
                sys.exit(f'\nError: bwa index failed to run on sample {s}:\n{p.stderr}')
        log(f'\r{i+1} / {len(samples)}', end='')
    log('\n')


def index_exists(in_dir, sample):
    index_1 = in_dir / (sample + '_full.fasta.amb')
    index_2 = in_dir / (sample + '_full.fasta.ann')
    index_3 = in_dir / (sample + '_full.fasta.bwt')
    index_4 = in_dir / (sample + '_full.fasta.pac')
    index_5 = in_dir / (sample + '_full.fasta.sa')
    return (index_1.is_file() and index_1.stat().st_size > 0 and
            index_2.is_file() and index_2.stat().st_size > 0 and
            index_3.is_file() and index_3.stat().st_size > 0 and
            index_4.is_file() and index_4.stat().st_size > 0 and
            index_5.is_file() and index_5.stat().st_size > 0)


def align_all_samples(in_dir, samples):
    section_header('Aligning pairwise combinations')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')
    for a in samples:
        for b in samples:
            align_sample_pair(in_dir, a, b)


def align_sample_pair(in_dir, a, b):
    pieces_fasta = in_dir / (a + '_pieces.fasta')
    sequence_fasta = in_dir / (b + '_full.fasta')
    assert index_exists(in_dir, b)
    log(f'{a} {b}')
