"""
Copyright 2022 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/XXXXXXXXX

This file is part of XXXXXXXXX. XXXXXXXXX is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. XXXXXXXXX is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with XXXXXXXXX.
If not, see <https://www.gnu.org/licenses/>.
"""

import pandas as pd
from plotnine import ggplot, aes, geom_segment, geom_line, geom_hline, geom_vline, labs, theme_bw, \
    scale_x_continuous, scale_x_sqrt, scale_y_continuous, scale_y_sqrt, scale_color_manual, \
    element_blank, theme

from .distance import get_distance


def show_plots(sample_name_a, sample_name_b, window_size, aligned_frac, masses, smoothed_masses,
               thresholds, painted_a, painted_b, sqrt_distance, sqrt_mass):
    distribution_plot(sample_name_a, sample_name_b, window_size, aligned_frac, masses,
                      smoothed_masses, thresholds, sqrt_distance, sqrt_mass)
    # contig_plot(sample_name_a, sample_name_b, aligned_frac, painted_a, window_size, thresholds)


def distribution_plot(sample_name_a, sample_name_b, window_size, aligned_frac, masses,
                      smoothed_masses, thresholds, sqrt_distance, sqrt_mass):
    title = f'{sample_name_a} vs {sample_name_b}, {window_size} bp windows, ' \
            f'{100.0 * aligned_frac:.1f}% aligned'

    mean = get_distance(masses, window_size, 'mean')
    # peak_mean = get_distance([m if low <= i <= high else 0.0 for i, m in enumerate(masses)],
    #                          window_size, 'mean')

    x_max = len(masses) / window_size
    y_max = 1.05 * max(max(masses), max(smoothed_masses))

    distances = [i / window_size for i in range(len(masses))]

    grouping = []
    low_2, low_1 = thresholds['low2'], thresholds['low1']
    high_1, high_2 = thresholds['high1'], thresholds['high2']
    for i in range(len(masses)):
        if low_2 is not None and i < low_2:
            grouping.append('low')
        elif low_1 is not None and i < low_1:
            grouping.append('lowish')
        elif high_2 is not None and i > high_2:
            grouping.append('high')
        elif high_1 is not None and i > high_1:
            grouping.append('highish')
        else:
            grouping.append('central')

    df = pd.DataFrame(list(zip(distances, masses, smoothed_masses, grouping)),
                      columns=['distance', 'mass', 'smoothed_mass', 'grouping'])

    g = (ggplot(df) +
         geom_segment(aes(x='distance', xend='distance', y=0, yend='mass', colour='grouping'),
                      size=1) +
         scale_color_manual({'low': '#d6d6d6', 'lowish': '#d7cfff',
                             'high': '#d6d6d6', 'highish': '#d7cfff',
                             'central': '#7570b3'}, guide=False) +
         geom_line(aes(x='distance', y='smoothed_mass'), size=0.5) +
         # geom_vline(xintercept=mean, colour='#d95f02', linetype='dotted', size=0.5) +
         # geom_vline(xintercept=peak_mean, colour='#d95f02', linetype='dashed', size=0.5) +
         theme_bw() +
         labs(title=title, x='distance', y='probability mass'))

    if sqrt_distance:
        g += scale_x_sqrt(limits=(0, x_max))
    else:
        g += scale_x_continuous(limits=(0, x_max))

    if sqrt_mass:
        g += scale_y_sqrt(expand=(0, 0), limits=(0, y_max))
    else:
        g += scale_y_continuous(expand=(0, 0), limits=(0, y_max))

    g.draw(show=True)


def contig_plot(sample_name_a, sample_name_b, aligned_frac, painted, window_size, thresholds):
    title = f'{sample_name_a} vs {sample_name_b}, {window_size} bp windows, ' \
            f'{100.0 * aligned_frac:.1f}% aligned'

    boundaries = [0]
    x_max = 0
    for name, contig in painted.contigs.items():
        x_max += contig.length
        boundaries.append(x_max)
    y_max = painted.get_max_differences()

    g = (ggplot() +
         theme_bw() +
         theme(panel_grid_major_x=element_blank(), panel_grid_minor_x=element_blank()) +
         geom_hline(yintercept=low, colour='#d95f02', linetype='dotted') +
         geom_hline(yintercept=high, colour='#d95f02', linetype='dotted') +
         scale_y_continuous(expand=(0, 0), limits=(0, y_max)) +
         scale_x_continuous(expand=(0, 0), limits=(0, x_max)) +
         labs(title=title, x='contig position', y='distance'))

    for b in boundaries:
        g += geom_vline(xintercept=b, colour='#aaaaaa', size=0.5)

    offset = 0
    for name, contig in painted.contigs.items():
        positions = [offset + d[0] for d in contig.window_differences]  # TEMP
        differences = [d[2] for d in contig.window_differences]  # TEMP
        df = pd.DataFrame(list(zip(positions, differences)), columns=['positions', 'differences'])
        g += geom_line(data=df, mapping=aes(x='positions', y='differences'), size=0.5)
        offset += contig.length

    g.draw(show=True)
