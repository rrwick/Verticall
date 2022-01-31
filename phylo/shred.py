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
import sys

from .log import log, section_header, explanation
from .misc import get_generator, get_sequence_file_type


def shred(args):
    welcome_message()
    assemblies = find_assemblies(args.in_dir, args.recursive)
    check_assemblies_and_out_dir(assemblies, args.out_dir)
    shred_all_assemblies(assemblies, args.size, args.overlap, args.out_dir)
    finished_message()


def welcome_message():
    section_header('Starting XXXXXXXXX shred')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')


def finished_message():
    section_header('Finished!')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')


def find_assemblies(in_dir, recursive):
    if recursive:
        files = [f for f in in_dir.glob('**/*') if f.is_file()]
    else:
        files = [f for f in in_dir.glob('*') if f.is_file()]
    assemblies = []
    for f in files:
        file_type = get_sequence_file_type(f)
        if file_type == 'FASTA' or file_type == 'GFA':
            assemblies.append(f)
    log(f'Found {len(assemblies):,} assembly files in {in_dir.resolve()}')
    log()
    return assemblies


def shred_all_assemblies(assemblies, size, overlap, out_dir):
    section_header('Shredding assemblies')
    explanation('Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor '
                'incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis '
                'nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.')
    for a in assemblies:
        shred_one_assembly(a, size, overlap, out_dir)
    log()


def shred_one_assembly(a, size, overlap, out_dir):
    assembly_name = get_name_from_filename(a)
    log(f'{assembly_name}: ')
    full_output, piece_output = get_output_filenames(assembly_name, out_dir)
    assembly_size, piece_count = 0, 0
    contig_names = set()
    with open(full_output, 'wt') as full_out, open(piece_output, 'wt') as piece_out:
        for contig_name, seq in get_generator(a)(a):
            if contig_name in contig_names:
                sys.exit(f'\nError: duplicate contig name in {a} ({contig_name})')
            contig_names.add(contig_name)
            full_out.write(f'>{contig_name}\n{seq}\n')
            assembly_size += len(seq)
            pieces = shred_sequence(seq, size, overlap)
            piece_count += len(pieces)
            for i, piece_seq in enumerate(pieces):
                piece_name = f'{assembly_name}_{contig_name}_{i:09}'
                piece_out.write(f'>{piece_name}\n{piece_seq}\n')
    log(f'  {full_output}:   {assembly_size:,} bp')
    log(f'  {piece_output}: {piece_count:,} pieces')
    log()


def shred_sequence(seq, size, overlap):
    """
    This function takes a sequence and returns a list of sequence pieces with a length defined by
    'size' and the length of shared sequences between adjacent pieces is defined by 'overlap'. In
    cases where these pieces don't entirely cover the sequence, the pieces will be centred in the
    sequence (i.e. the unused parts of the sequence will be at both the start and end).
    """
    pieces = []
    if len(seq) < size:
        return pieces

    step = size - overlap
    steps_in_seq = (len(seq) - size) // step
    used_len = size + (step * steps_in_seq)
    seq = trim_seq(seq, used_len)

    start, end = 0, size
    while end <= len(seq):
        piece = seq[start:end]
        assert len(piece) == size
        pieces.append(piece)
        start += step
        end += step

    assert len(pieces) == steps_in_seq + 1
    return pieces


def trim_seq(seq, used_len):
    trim_amount = len(seq) - used_len
    assert trim_amount >= 0  # can't make the sequence longer
    if trim_amount == 0:
        return seq
    start_trim = trim_amount // 2
    if trim_amount % 2 == 0:  # even trim - same amount from start and end
        end_trim = start_trim
    else:  # odd trim - one more from end
        end_trim = start_trim + 1
    assert start_trim + end_trim + used_len == len(seq)
    return seq[start_trim:-end_trim]


def check_assemblies_and_out_dir(assemblies, out_dir):
    if out_dir.is_dir():
        log('Output directory already exists:')
    else:
        log('Creating output directory:')
        out_dir.mkdir(parents=True, exist_ok=True)
    log(f'  {out_dir.resolve()}')
    log()

    assembly_names = set()
    for a in assemblies:
        name = get_name_from_filename(a)
        if name in assembly_names:
            sys.exit(f'\nError: duplicate assembly name in input directory ({name})')
        assembly_names.add(name)
        full_output, piece_output = get_output_filenames(name, out_dir)
        if full_output.is_file():
            sys.exit(f'\nError: output file already exists ({full_output})')
        if piece_output.is_file():
            sys.exit(f'\nError: output file already exists ({piece_output})')


def get_name_from_filename(f: pathlib.Path):
    """
    Takes an assembly filepath and returns the assembly name. First '.gz' is dropped (if present)
    and then the first file extension.
    """
    if str(f).endswith('.gz'):
        return pathlib.Path(str(f)[:-3]).stem
    else:
        return f.stem


def get_output_filenames(name, out_dir):
    full_output = out_dir / (name + '_full.fasta')
    piece_output = out_dir / (name + '_pieces.fasta')
    return full_output, piece_output
