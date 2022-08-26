#!/usr/bin/env python3
"""
This script takes in two whole-genome pseudo-alignments (before and after masking), and it produces
an SVG image showing masked and unmasked regions of each sequence. It uses three colours:
* Unmasked
* Missing (regions that are missing in the unmasked sequence)
* Masked (regions that were masked by Gubbins/ClonalFrameML/Verticall)

It has three required arguments:
* --unmasked: the full unmasked alignment
* --masking: the masking information, which can be:
  * Gubbins' recombination_predictions.gff file
  * ClonalFrameML's importation_status.txt and labelled_tree.newick files (comma-delimited)
  * Verticall's TSV file
  * 'None' (will only visualise the unmasked alignment)

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
import newick
import pathlib
import re
import svgwrite
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Draw Gubbins masking')

    required_args = parser.add_argument_group('Required')
    required_args.add_argument('--unmasked', type=str, required=True,
                             help='Unmasked alignment')
    required_args.add_argument('--masking', type=str, required=True,
                               help='Masking info')
    required_args.add_argument('--image', type=str, required=True,
                               help='Output filename of SVG illustration')

    order_args = parser.add_argument_group('Sample order')
    order_args.add_argument('--tree', type=str,
                            help='Tree in newick format to specify order of samples')
    order_args.add_argument('--reverse', action="store_true",
                            help='Reverse the order of samples')

    spacing_args = parser.add_argument_group('Spacing')
    spacing_args.add_argument('--line_width', type=float, default=8,
                              help='Width of lines')
    spacing_args.add_argument('--gap', type=float, default=1,
                              help='Gap between lines')

    colour_args = parser.add_argument_group('Colours')
    colour_args.add_argument('--unmasked_colour', type=str, default='#cccccc',
                             help='Hex colour for unmasked regions')
    colour_args.add_argument('--missing_colour', type=str, default='#eeeeee',
                             help='Hex colour for missing regions')
    colour_args.add_argument('--masked_colour', type=str, default='#9b0000',
                             help='Hex colour for masked regions')
    colour_args.add_argument('--ignore_missing_size', type=int, default=10,
                             help="Don't draw missing regions this size or smaller")


    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    alignment_names, alignment_seqs, alignment_length = load_alignment(args.unmasked)
    names = get_names(alignment_names, args.tree, args.reverse)
    masked_regions, mask_type = get_masked_regions(args.masking)

    image = svgwrite.Drawing(args.image, profile='full')
    y_pos = args.line_width + args.gap

    for name in names:
        print(name)
        image.add(image.text(name, insert=(97, y_pos+(args.line_width/3.0)),
                             style='text-anchor:end', font_size=f'{args.line_width}px'))
        image.add(image.line((100, y_pos), (500, y_pos), stroke=args.unmasked_colour,
                             stroke_width=args.line_width))
        missing_regions = get_missing_regions(alignment_seqs[name], args.ignore_missing_size)
        if mask_type == 'Verticall':
            verticall_horizontal, verticall_unaligned = masked_regions
            horizontal_regions = verticall_horizontal[name]
            missing_regions += verticall_unaligned[name]
        else:
            horizontal_regions = masked_regions[name]

        for start, end in missing_regions:
            image.add(image.line((100 + 400 * start / alignment_length, y_pos),
                                 (100 + 400 * end / alignment_length, y_pos),
                                 stroke=args.missing_colour, stroke_width=args.line_width))
        for start, end in horizontal_regions:
            print(f'  {start}-{end}')
            image.add(image.line((100 + 400 * start / alignment_length, y_pos),
                                 (100 + 400 * end / alignment_length, y_pos),
                                 stroke=args.masked_colour, stroke_width=args.line_width))
        y_pos += args.line_width
        y_pos += args.gap
    image.save()


def get_names(alignment_names, tree_filename, reverse):
    if tree_filename is not None:
        tree = newick.read(tree_filename)[0]
        names = tree.get_leaf_names()
        alignment_names = set(alignment_names)
        for name in names:
            if name not in alignment_names:
                sys.exit(f'Error: {name} in tree but not in alignment')
    else:
        names = alignment_names
    if reverse:
        names = names[::-1]
    return names


def get_missing_regions(unmasked_seq, ignore_size):
    missing_regions = []
    start = None
    for i, unmasked_base in enumerate(unmasked_seq):
        if not is_canonical(unmasked_base):
            if start is None:
                start = i
        else:
            if start is not None:
                if i - start > ignore_size:
                    missing_regions.append((start, i))
                start = None
    if start is not None:
        i += 1
        if i - start > ignore_size:
            missing_regions.append((start, i))
    return missing_regions


def get_masked_regions(masking_file):
    if masking_file == 'None':
        return collections.defaultdict(list), 'None'
    elif masking_file.endswith('.tsv'):
        return get_masked_regions_verticall(masking_file), 'Verticall'
    elif masking_file.endswith('.gff'):
        return get_masked_regions_gubbins(masking_file), 'Gubbins'
    elif ',' in masking_file:
        return get_masked_regions_clonalframeml(masking_file), 'ClonalFrameML'
    else:
        sys.exit('Error: could not decipher masking files')


def get_masked_regions_verticall(tsv_filename):
    print(f'Loading Verticall masked regions from {tsv_filename}')
    masked_regions_horizontal = collections.defaultdict(list)
    masked_regions_unaligned = collections.defaultdict(list)
    name_col, h_col, u_col = None, None, None
    with open(tsv_filename, 'rt') as tsv_file:
        for i, line in enumerate(tsv_file):
            parts = line.strip('\n').split('\t')
            if i == 0:
                name_col = get_column_index(parts, 'assembly_b')
                h_col = get_column_index(parts, 'assembly_a_horizontal_regions')
                u_col = get_column_index(parts, 'assembly_a_unaligned_regions')
                continue
            name = parts[name_col]
            if name in masked_regions_horizontal:  # seen this one already
                continue
            horizontal_regions = [parse_region(x) for x in parts[h_col].split(',') if x]
            unaligned_regions = [parse_region(x) for x in parts[u_col].split(',') if x]
            masked_regions_horizontal[name] = horizontal_regions
            masked_regions_unaligned[name] = unaligned_regions
    return masked_regions_horizontal, masked_regions_unaligned


def get_column_index(header_parts, column):
    for i, header in enumerate(header_parts):
        if column == header:
            return i
    sys.exit(f'Error: could not find {column} column')


def parse_region(region_str):
    name, region = region_str.split(':')
    start, end = region.split('-')
    return int(start), int(end)


def get_masked_regions_gubbins(gff_filename):
    """
    This function uses some code from the mask_gubbins_aln.py script in Gubbins:
    https://github.com/nickjcroucher/gubbins/blob/master/python/scripts/mask_gubbins_aln.py
    """
    print(f'Loading Gubbins masked regions from {gff_filename}')
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
    for taxon, regions in masked_regions.items():
        masked_regions[taxon] = sorted(regions)
    return masked_regions


def get_masked_regions_clonalframeml(clonalframeml_files):
    clonalframeml_files = clonalframeml_files.split(',')
    assert len(clonalframeml_files) == 2
    if clonalframeml_files[0].endswith('.txt') and clonalframeml_files[1].endswith('.newick'):
        imports_filename, tree_filename = clonalframeml_files
    elif clonalframeml_files[0].endswith('.newick') and clonalframeml_files[1].endswith('.txt'):
        tree_filename, imports_filename = clonalframeml_files
    else:
        sys.exit('Error: could not decipher masking files')
    print(f'Loading ClonalFrameML masked regions from {tree_filename} and {imports_filename}')
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


def is_canonical(b):
    return b == 'A' or b == 'C' or b == 'G' or b == 'T'


def load_alignment(alignment_filename):
    names = []
    seqs = {}
    alignment_length = None
    for name, seq in iterate_fasta(alignment_filename):
        names.append(name)
        seqs[name] = seq
        if alignment_length is None:
            alignment_length = len(seq)
        else:
            assert alignment_length == len(seq)
    return names, seqs, alignment_length


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


# Unit tests for Pytest
# =====================

def test_get_missing_regions():
    assert get_missing_regions('ACGTACGATC', 0) == []
    assert get_missing_regions('ANGTACGATC', 0) == [(1, 2)]
    assert get_missing_regions('ANNNACGANN', 0) == [(1, 4), (8, 10)]
    assert get_missing_regions('NNNNNNNNNN', 0) == [(0, 10)]

    assert get_missing_regions('ACGTACGATC', 1) == []
    assert get_missing_regions('ANGTACGATC', 1) == []
    assert get_missing_regions('ANNNACGANN', 1) == [(1, 4), (8, 10)]
    assert get_missing_regions('NNNNNNNNNN', 1) == [(0, 10)]

    assert get_missing_regions('ACGTACGATC', 2) == []
    assert get_missing_regions('ANGTACGATC', 2) == []
    assert get_missing_regions('ANNNACGANN', 2) == [(1, 4)]
    assert get_missing_regions('NNNNNNNNNN', 2) == [(0, 10)]

    assert get_missing_regions('ACGTACGATC', 3) == []
    assert get_missing_regions('ANGTACGATC', 3) == []
    assert get_missing_regions('ANNNACGANN', 3) == []
    assert get_missing_regions('NNNNNNNNNN', 3) == [(0, 10)]

    assert get_missing_regions('ACGTACGATC', 9) == []
    assert get_missing_regions('ANGTACGATC', 9) == []
    assert get_missing_regions('ANNNACGANN', 9) == []
    assert get_missing_regions('NNNNNNNNNN', 9) == [(0, 10)]

    assert get_missing_regions('ACGTACGATC', 10) == []
    assert get_missing_regions('ANGTACGATC', 10) == []
    assert get_missing_regions('ANNNACGANN', 10) == []
    assert get_missing_regions('NNNNNNNNNN', 10) == []


def test_get_masked_regions():
    assert get_masked_regions('ACGTACGATC', 'ACGTACGATC', 0) == []
    assert get_masked_regions('ANNNACGANN', 'ANNNACGANN', 0) == []
    assert get_masked_regions('NNNNNNNNNN', 'NNNNNNNNNN', 0) == []

    assert get_masked_regions('ACGTACGATC', 'ANGTACGATC', 0) == [(1, 2)]
    assert get_masked_regions('ACGTACGATC', 'ANNNACGANN', 0) == [(1, 4), (8, 10)]
    assert get_masked_regions('ACGTACGATC', 'NNNNNNNNNN', 0) == [(0, 10)]

    assert get_masked_regions('ACGTACGATC', 'ANGTACGATC', 1) == []
    assert get_masked_regions('ACGTACGATC', 'ANNNACGANN', 1) == [(1, 4), (8, 10)]
    assert get_masked_regions('ACGTACGATC', 'NNNNNNNNNN', 1) == [(0, 10)]

    assert get_masked_regions('ACGTACGATC', 'ANGTACGATC', 2) == []
    assert get_masked_regions('ACGTACGATC', 'ANNNACGANN', 2) == [(1, 4)]
    assert get_masked_regions('ACGTACGATC', 'NNNNNNNNNN', 2) == [(0, 10)]

    assert get_masked_regions('ACGTACGATC', 'ANGTACGATC', 3) == []
    assert get_masked_regions('ACGTACGATC', 'ANNNACGANN', 3) == []
    assert get_masked_regions('ACGTACGATC', 'NNNNNNNNNN', 3) == [(0, 10)]

    assert get_masked_regions('ACGTACGATC', 'ANGTACGATC', 9) == []
    assert get_masked_regions('ACGTACGATC', 'ANNNACGANN', 9) == []
    assert get_masked_regions('ACGTACGATC', 'NNNNNNNNNN', 9) == [(0, 10)]

    assert get_masked_regions('ACGTACGATC', 'ANGTACGATC', 10) == []
    assert get_masked_regions('ACGTACGATC', 'ANNNACGANN', 10) == []
    assert get_masked_regions('ACGTACGATC', 'NNNNNNNNNN', 10) == []
