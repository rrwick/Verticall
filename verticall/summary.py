#!/usr/bin/env python3
"""
This module contains code for the 'verticall summary' subcommand.

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

import matplotlib.pyplot as plt
import pandas as pd
from plotnine import ggplot, aes, geom_area, geom_vline, labs, theme_bw, scale_x_continuous, \
    scale_y_continuous, scale_fill_manual, element_blank, theme
import sys

from .matrix import get_column_index
from .misc import iterate_fasta, get_open_func
from .tsv import split_region_str


def summary(args):
    contig_lengths, sample_name = get_contig_lengths(args.assembly)
    data = load_data(args.in_file, sample_name)
    summarised_data = summarise_data(data, contig_lengths, args.all)
    if args.plot:
        plot = summary_plot(sample_name, summarised_data, contig_lengths, args.vertical_colour,
                            args.horizontal_colour, args.unaligned_colour)
        plt.show()
    else:
        print('contig', 'position', 'vertical', 'horizontal', 'unaligned')
        for contig, position, vertical, horizontal, unaligned in summarised_data:
            print(f'{contig}\t{position}\t{vertical}\t{horizontal}\t{unaligned}')


def get_contig_lengths(assembly_filename):
    extension_len = None
    if str(assembly_filename).endswith('.fasta'):
        extension_len = 6
    elif str(assembly_filename).endswith('.fasta.gz'):
        extension_len = 9
    elif str(assembly_filename).endswith('.fna'):
        extension_len = 4
    elif str(assembly_filename).endswith('.fna.gz'):
        extension_len = 7
    elif str(assembly_filename).endswith('.fa'):
        extension_len = 3
    elif str(assembly_filename).endswith('.fa.gz'):
        extension_len = 6
    if extension_len is None:
        sys.exit(f'Error: {assembly_filename} does not end in a FASTA file extension')
    sample_name = assembly_filename.name[:-extension_len]
    contig_lengths = {}
    for name, seq in iterate_fasta(assembly_filename):
        contig_lengths[name] = len(seq)
    return contig_lengths, sample_name


def load_data(filename, sample_name):
    data = []
    v_column, h_column, u_column = None, None, None
    with get_open_func(filename)(filename, 'rt') as pairwise_file:
        for i, line in enumerate(pairwise_file):
            parts = line.strip('\n').split('\t')
            if i == 0:
                v_column = get_column_index(parts, 'assembly_a_vertical_regions', filename)
                h_column = get_column_index(parts, 'assembly_a_horizontal_regions', filename)
                u_column = get_column_index(parts, 'assembly_a_unaligned_regions', filename)
            elif parts[0] == sample_name:
                vertical_regions = parts[v_column].split(',') if parts[v_column] else []
                horizontal_regions = parts[h_column].split(',') if parts[h_column] else []
                unaligned_regions = parts[u_column].split(',') if parts[u_column] else []
                data.append((vertical_regions, horizontal_regions, unaligned_regions))
    return data


def summarise_data(data, contig_lengths, output_all):
    vertical_counts = {name: [0] * length for name, length in contig_lengths.items()}
    horizontal_counts = {name: [0] * length for name, length in contig_lengths.items()}
    unaligned_counts = {name: [0] * length for name, length in contig_lengths.items()}
    for vertical_regions, horizontal_regions, unaligned_regions in data:
        for region in vertical_regions:
            name, start, end = split_region_str(region)
            for i in range(start, end):
                vertical_counts[name][i] += 1
        for region in horizontal_regions:
            name, start, end = split_region_str(region)
            for i in range(start, end):
                horizontal_counts[name][i] += 1
        for region in unaligned_regions:
            name, start, end = split_region_str(region)
            for i in range(start, end):
                unaligned_counts[name][i] += 1
    summarised_data = []
    for name, length in contig_lengths.items():
        for i in range(length):
            vertical_count = vertical_counts[name][i]
            horizontal_count = horizontal_counts[name][i]
            unaligned_count = unaligned_counts[name][i]
            if i == 0 or i == length-1 or output_all:
                summarised_data.append((name, i, vertical_count, horizontal_count, unaligned_count))
            else:
                prev_counts = (vertical_counts[name][i-1], horizontal_counts[name][i-1],
                               unaligned_counts[name][i-1])
                this_counts = (vertical_count, horizontal_count, unaligned_count)
                next_counts = (vertical_counts[name][i+1], horizontal_counts[name][i+1],
                               unaligned_counts[name][i+1])
                if prev_counts != this_counts or next_counts != this_counts:
                    summarised_data.append((name, i, vertical_count, horizontal_count,
                                            unaligned_count))
    return summarised_data


def summary_plot(sample_name, summarised_data, contig_lengths, vertical_colour, horizontal_colour,
                 unaligned_colour):
    title = f'{sample_name} painting summary'

    boundaries = [0]
    x_max = 0
    for name, length in contig_lengths.items():
        x_max += length
        boundaries.append(x_max)
    y_max = 0
    for _, _, vertical, horizontal, unaligned in summarised_data:
        y_max = max(y_max, vertical + horizontal + unaligned)

    g = (ggplot() +
         theme_bw() + theme(figure_size=(15, 5)) +
         theme(panel_grid_major_x=element_blank(), panel_grid_minor_x=element_blank()) +
         scale_x_continuous(expand=(0, 0), limits=(0, x_max)) +
         scale_y_continuous(expand=(0, 0), limits=(0, y_max)) +
         labs(title=title, x='contig position', y='count'))

    df = pd.DataFrame(summarised_data,
                      columns=['contig', 'pos', 'vertical', 'horizontal', 'unaligned'])

    # Make the data tidy.
    df = pd.melt(df, id_vars=['contig', 'pos'],
                 value_vars=['vertical', 'horizontal', 'unaligned'],
                 var_name='classification', value_name='count')
    df['classification'] = df['classification'].astype('category')
    df['classification'] = \
        df['classification'].cat.reorder_categories(['unaligned', 'horizontal', 'vertical'])

    offset = 0
    for name, length in contig_lengths.items():
        contig_df = df[df['contig'] == name].copy()
        contig_df['offset_pos'] = contig_df['pos'] + offset
        g += geom_area(data=contig_df,
                       mapping=aes(x='offset_pos', y='count', fill='classification'))
        offset += length

    g += scale_fill_manual({'vertical': vertical_colour, 'horizontal': horizontal_colour,
                            'unaligned': unaligned_colour}, guide=False)

    for b in boundaries:
        g += geom_vline(xintercept=b, colour='#000000', size=0.5)

    return g.draw()
