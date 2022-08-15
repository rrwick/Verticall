#!/usr/bin/env python3
"""
This script produces a masked ClonalFrameML alignment. It takes as input an unmasked whole-genome
pseudo-alignment, ClonalFrameML's importation_status.txt file (which describes the recombination
regions) and ClonalFrameML's tree. It outputs (to stdout) a masked version of the alignment, with
Ns in place of masked sequences.

Example usage:
  mask_clonalframeml_alignment.py \
    alignment.fasta \
    clonalframeml.importation_status.txt \
    clonalframeml.labelled_tree.newick > masked_alignment.fasta

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
import newick
import gzip
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Mask ClonalFrameML alignment')

    parser.add_argument('unmasked', type=str,
                        help='Filename of unmasked alignment in FASTA format')
    parser.add_argument('imports', type=str,
                        help="Filename of ClonalFrameML's importation_status.txt file")
    parser.add_argument('tree', type=str,
                        help="Filename of ClonalFrameML's labelled_tree.newick file")

    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    masked_regions = load_clonalframeml_regions(args.imports, args.tree)
    for name, seq in iterate_fasta(args.unmasked):
        seq = mask_sequence(seq, masked_regions[name])
        print(f'>{name}\n{seq}')


def load_clonalframeml_regions(imports_filename, tree_filename):
    imports = load_imports(imports_filename)
    tree = newick.read(tree_filename)[0]
    masked_regions = collections.defaultdict(list)
    for node, start, end in imports:
        tips = tree.get_node(node).get_leaf_names()
        for tip in tips:
            masked_regions[tip].append((start, end))
    return masked_regions


def load_imports(imports_filename):
    imports = []
    with open(imports_filename, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts[1] == 'Beg':  # header line
                continue
            node, start, end = parts[0], int(parts[1]), int(parts[2])
            start -= 1  # change from 1-based inclusive to 0-based exclusive
            imports.append((node, start, end))
    return imports


def mask_sequence(seq, masked_regions):
    seq = [base for base in seq]
    for start, end in masked_regions:
        assert start >= 0 and end <= len(seq)
        for i in range(start, end):
            if seq[i] != '-':
                seq[i] = 'N'
    return ''.join(seq)


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def iterate_fasta(filename, preserve_case=False):
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
                if preserve_case:
                    sequence.append(line)
                else:
                    sequence.append(line.upper())
        if name:
            yield name.split()[0], ''.join(sequence)


if __name__ == '__main__':
    main()
