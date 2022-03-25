#!/usr/bin/env python3
"""
This is a script which takes two arguments:
* the pairwise.tsv file made by Verticall
* a sample name

And it outputs (to stdout) a table describing (vertical, horizontal or unaligned) each position
of the assembly.

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
import matplotlib.pyplot as plt
import pandas as pd
from plotnine import ggplot, aes, geom_area, geom_vline, labs, theme_bw, scale_x_continuous, \
    scale_y_continuous, scale_fill_manual, element_blank, theme

VERTICAL_COLOUR = '#4859a0'
HORIZONTAL_COLOUR = '#c47e7e'
UNALIGNED_COLOUR = '#eeeeee'


def get_arguments():
    parser = argparse.ArgumentParser(description='Summarise positions of an assembly')

    parser.add_argument('pairwise', type=str,
                        help='Vertical\'s pairwise.tsv file')
    parser.add_argument('name', type=str,
                        help='Sample name for the assembly to summarise')

    parser.add_argument('--all', action='store_true',
                        help='Output one line for all assembly positions (default: omit redundant '
                             'adjacent lines)')
    parser.add_argument('--plot', action='store_true',
                        help='Instead of outputting a table, display an interactive plot')
    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    data = load_data(args.pairwise, args.name)
    contig_lengths = get_contig_lengths(data)
    summarised_data = summarise_data(data, contig_lengths, args.all)
    if args.plot:
        plot = summary_plot(args.name, summarised_data, contig_lengths)
        plt.show()
    else:
        print('contig', 'position', 'vertical', 'horizontal', 'unaligned')
        for contig, position, vertical, horizontal, unaligned in summarised_data:
            print(f'{contig}\t{position}\t{vertical}\t{horizontal}\t{unaligned}')


def load_data(pairwise_filename, sample_name):
    data = []
    with open(pairwise_filename, 'rt') as pairwise_file:
        for i, line in enumerate(pairwise_file):
            if i == 0:
                continue  # skip header line
            parts = line.strip().split('\t')
            assert len(parts) == 27
            if parts[0] == sample_name:
                vertical_regions = parts[21].split(',') if parts[21] else []
                horizontal_regions = parts[22].split(',') if parts[22] else []
                unaligned_regions = parts[23].split(',') if parts[23] else []
                data.append((vertical_regions, horizontal_regions, unaligned_regions))
    return data


def get_contig_lengths(data):
    contig_lengths = collections.defaultdict(int)
    for vertical_regions, horizontal_regions, unaligned_regions in data:
        all_regions = vertical_regions + horizontal_regions + unaligned_regions
        for region in all_regions:
            name, start, end = split_region_str(region)
            contig_lengths[name] = max(contig_lengths[name], int(end))
    return contig_lengths


def split_region_str(region):
    name, contig_range = region.split(':')
    start, end = contig_range.split('-')
    return name, int(start), int(end)


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


def summary_plot(sample_name, summarised_data, contig_lengths):
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
         theme_bw() +
         theme(panel_grid_major_x=element_blank(), panel_grid_minor_x=element_blank()) +
         scale_x_continuous(expand=(0, 0), limits=(0, x_max)) +
         scale_y_continuous(expand=(0, 0), limits=(0, y_max)) +
         labs(title=title, x='contig position', y='counts'))

    df = pd.DataFrame(summarised_data,
                      columns=['contig', 'pos', 'vertical', 'horizontal', 'unaligned'])

    # Make the data tidy.
    df = pd.melt(df, id_vars=['contig', 'pos'],
                 value_vars=['vertical', 'horizontal', 'unaligned'],
                 var_name='classification', value_name='count')

    offset = 0
    for name, length in contig_lengths.items():
        contig_df = df[df['contig'] == name]
        contig_df['pos'] += offset
        g += geom_area(data=contig_df, mapping=aes(x='pos', y='count', fill='classification'))
        offset += length

    g += scale_fill_manual({'vertical': VERTICAL_COLOUR, 'horizontal': HORIZONTAL_COLOUR,
                            'unaligned': UNALIGNED_COLOUR}, guide=False)

    for b in boundaries:
        g += geom_vline(xintercept=b, colour='#000000', size=0.5)

    return g.draw()


if __name__ == '__main__':
    main()
