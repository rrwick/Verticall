#!/usr/bin/env python3
"""
This script takes a Gubbins GFF file as input and produces an SVG image in the same style as
Verticall mask. It can be useful for comparing which regions of an alignment were masked by Gubbins
vs Verticall.

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
import pathlib
import svgwrite
import re
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Draw Gubbins masking')

    parser.add_argument('alignment', type=pathlib.Path,
                        help='Input alignment')
    parser.add_argument('gff', type=pathlib.Path,
                        help='Input filename of Gubbins GFF file')
    parser.add_argument('image', type=pathlib.Path,
                        help='Output filename of SVG illustration of masked regions')

    parser.add_argument('--vertical_colour', type=str, default='#4859a0',
                        help='Hex colour for vertical inheritance')
    parser.add_argument('--horizontal_colour', type=str, default='#c47e7e',
                        help='Hex colour for horizontal inheritance')

    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    names, alignment_length = load_alignment(args.alignment)
    masked_regions = load_gubbins_masked_regions(args.gff)
    image = svgwrite.Drawing(args.image, profile='full')
    y_pos = 12
    for sample_name in names:
        print(f'{sample_name}:')
        image.add(image.text(sample_name, insert=(97, y_pos+4), style='text-anchor:end',
                             font_size='12px'))
        image.add(image.line((100, y_pos), (500, y_pos), stroke=args.vertical_colour,
                             stroke_width=9))
        sample_name = sample_name.replace('#', '_')
        if sample_name in masked_regions:
            for start, end in masked_regions[sample_name]:
                print(f'  {start}-{end}')
                image.add(image.line((100 + 400 * start / alignment_length, y_pos),
                                     (100 + 400 * end / alignment_length, y_pos),
                                     stroke=args.horizontal_colour, stroke_width=9))
        y_pos += 12
    image.save()


def load_alignment(alignment_filename):
    names = []
    alignment_length = None
    for name, seq in iterate_fasta(alignment_filename):
        names.append(name)
        if alignment_length is None:
            alignment_length = len(seq)
        else:
            assert alignment_length == len(seq)
    return names, alignment_length


def load_gubbins_masked_regions(gff_filename):
    """
    This function uses some code from the mask_gubbins_aln.py script in Gubbins:
    https://github.com/nickjcroucher/gubbins/blob/master/python/scripts/mask_gubbins_aln.py
    """
    masked_regions = collections.defaultdict(list)
    taxon_pattern = re.compile('taxa="([^"]*)"')
    with get_open_func(gff_filename)(gff_filename, 'r') as gff_file:
        for line in gff_file.readlines():
            if not line.startswith('##'):
                info = line.rstrip().split('\t')
                start = int(info[3]) - 1
                end = int(info[4]) - 1
                taxa = set(taxon_pattern.search(info[8]).group(1).split())
                for taxon in taxa:
                    masked_regions[taxon].append((start, end))
    return {taxon: sorted(regions) for taxon, regions in masked_regions.items()}


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
