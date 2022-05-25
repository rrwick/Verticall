#!/usr/bin/env python3
"""
This is a script which takes one argument: a FASTA whole-genome pseudo-alignment.

And it outputs (to stdout) a single sequence that contains the consensus base at each position.

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
import collections
import gzip
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Produce consensus sequence from alignment')

    parser.add_argument('alignment', type=str,
                        help='Filename of whole-genome pseudo-alignment')
    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    alignment = {name: seq for name, seq in iterate_fasta(args.alignment)}
    seq = get_consensus_sequence(alignment)
    print(f'>consensus\n{seq}')


def get_consensus_sequence(sequences):
    """
    Returns an alignment where any columns that lack variation are removed.
    """
    alignment_length = get_alignment_length(sequences)
    unambiguous_bases = {'A', 'C', 'G', 'T'}
    consensus_seq = []
    for i in range(alignment_length):
        bases_at_pos = [seq[i] for seq in sequences.values() if seq[i] in unambiguous_bases]
        consensus_seq.append(most_common(bases_at_pos))
        if i % 1000 == 0:
            print(f'\r{i:,} / {alignment_length:,}', end='', flush=True, file=sys.stderr)
    print(f'\r{alignment_length:,} / {alignment_length:,}', file=sys.stderr)
    return ''.join(consensus_seq)


def most_common(lst):
    """
    https://stackoverflow.com/questions/1518522/find-the-most-common-element-in-a-list
    """
    data = collections.Counter(lst)
    return max(lst, key=data.get)


def get_alignment_length(sequences):
    alignment_lengths = {len(seq) for seq in sequences.values()}
    assert len(alignment_lengths) == 1
    return list(alignment_lengths)[0]


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


def iterate_fasta(filename):
    """
    Takes a FASTA file as input and yields the contents as (name, seq) tuples.
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
                    yield name.split()[0], ''.join(sequence)
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line.upper())
        if name:
            yield name.split()[0], ''.join(sequence)


if __name__ == '__main__':
    main()
