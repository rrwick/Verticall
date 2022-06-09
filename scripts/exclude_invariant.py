#!/usr/bin/env python3
"""
This is a script which takes one argument: a FASTA whole-genome pseudo-alignment.

And it outputs (to stdout) an alignment with all invariant sites removed.

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

import argparse
import gzip
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Exclude invariant sites in alignment')

    parser.add_argument('alignment', type=str,
                        help='Filename of whole-genome pseudo-alignment')
    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    alignment = {name: seq for name, seq in iterate_fasta(args.alignment, preserve_case=True)}
    alignment = drop_invariant_positions(alignment)
    for name, seq in alignment.items():
        print(f'>{name}\n{seq}')
    print(file=sys.stderr)


def drop_invariant_positions(sequences):
    """
    Returns an alignment where any columns that lack variation are removed.
    """
    alignment_length = get_alignment_length(sequences)
    positions_to_remove = set()
    a, c, g, t, n = 0, 0, 0, 0, 0
    print(f'0 / {alignment_length:,}', file=sys.stderr, end='', flush=True)
    for i in range(alignment_length):
        if (i+1) % 1000 == 0:
            print(f'\r{i+1:,} / {alignment_length:,}', file=sys.stderr, end='', flush=True)
        bases_at_pos = {seq[i] for seq in sequences.values()}
        real_base_count = count_real_bases(bases_at_pos)
        if real_base_count == 1:
            positions_to_remove.add(i)
            if 'A' in bases_at_pos or 'a' in bases_at_pos:
                a += 1
            elif 'C' in bases_at_pos or 'c' in bases_at_pos:
                c += 1
            elif 'G' in bases_at_pos or 'g' in bases_at_pos:
                g += 1
            elif 'T' in bases_at_pos or 't' in bases_at_pos:
                t += 1
            else:
                assert False
        elif real_base_count == 0:
            positions_to_remove.add(i)
            n += 1
    print(f'\r{alignment_length:,} / {alignment_length:,}', file=sys.stderr, end='\n\n', flush=True)
    assert a + c + g + t + n == len(positions_to_remove)
    if not positions_to_remove:
        print(f'no invariant positions removed from pseudo-alignment', file=sys.stderr)
    else:
        percentage = 100.0 * len(positions_to_remove)/alignment_length
        print(f'{len(positions_to_remove):,} invariant positions ({percentage:.3}%) removed from '
              f'pseudo-alignment:', file=sys.stderr)
        print(f'  {a:9,} × A', file=sys.stderr)
        print(f'  {c:9,} × C', file=sys.stderr)
        print(f'  {g:9,} × G', file=sys.stderr)
        print(f'  {t:9,} × T', file=sys.stderr)
        print(f'  {n:9,} × other', file=sys.stderr)
    return drop_positions(sequences, positions_to_remove)


def get_alignment_length(sequences):
    alignment_lengths = {len(seq) for seq in sequences.values()}
    assert len(alignment_lengths) == 1
    return list(alignment_lengths)[0]


def drop_positions(sequences, positions_to_remove):
    if len(positions_to_remove) == 0:
        return sequences
    new_sequences = {}
    for name, seq in sequences.items():
        new_seq = ''.join(b for i, b in enumerate(seq) if i not in positions_to_remove)
        new_sequences[name] = new_seq
    return new_sequences


def count_real_bases(base_set):
    count = 0
    if 'A' in base_set or 'a' in base_set:
        count += 1
    if 'C' in base_set or 'c' in base_set:
        count += 1
    if 'G' in base_set or 'g' in base_set:
        count += 1
    if 'T' in base_set or 't' in base_set:
        count += 1
    return count


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    https://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    unknown_file = open(str(filename), 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('\nError: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('\nError: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def iterate_fasta(filename, include_info=False, preserve_case=False):
    """
    Takes a FASTA file as input and yields the contents as (name, seq) tuples. If include_info is
    set, it will yield (name, info, seq) tuples, where info is whatever follows the name.
    """
    with get_open_func(filename)(filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    if include_info:
                        name_parts = name.split(maxsplit=1)
                        contig_name = name_parts[0]
                        info = '' if len(name_parts) == 1 else name_parts[1]
                        yield contig_name, info, ''.join(sequence)
                    else:
                        yield name.split()[0], ''.join(sequence)
                    sequence = []
                name = line[1:]
            else:
                if preserve_case:
                    sequence.append(line)
                else:
                    sequence.append(line.upper())
        if name:
            if include_info:
                name_parts = name.split(maxsplit=1)
                contig_name = name_parts[0]
                info = '' if len(name_parts) == 1 else name_parts[1]
                yield contig_name, info, ''.join(sequence)
            else:
                yield name.split()[0], ''.join(sequence)


if __name__ == '__main__':
    main()
