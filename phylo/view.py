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
from plotnine import ggplot, aes, geom_segment, geom_line, geom_vline, labs, theme_bw, \
    scale_x_continuous, scale_x_sqrt, scale_y_continuous, scale_y_sqrt, scale_color_manual
import sys

from .distance import get_distance


def show_plots(sample_name_a, sample_name_b, window_size, aligned_frac, masses, smoothed_masses,
               low, high, sqrt_x, sqrt_y):
    distribution_plot(sample_name_a, sample_name_b, window_size, aligned_frac, masses,
                      smoothed_masses, low, high, sqrt_x, sqrt_y)


def distribution_plot(sample_name_a, sample_name_b, window_size, aligned_frac, masses,
                      smoothed_masses, low, high, sqrt_x, sqrt_y):
    title = f'{sample_name_a} vs {sample_name_b}, {window_size} bp windows, ' \
            f'{100.0 * aligned_frac:.1f}% aligned'

    mean = get_distance(masses, window_size, 'mean')
    peak_mean = get_distance([m if low <= i <= high else 0.0 for i, m in enumerate(masses)],
                             window_size, 'mean')

    x_max = len(masses) / window_size
    y_max = 1.05 * max(max(masses), max(smoothed_masses))

    distances = [i / window_size for i in range(len(masses))]
    in_main_peak = [low <= i <= high for i in range(len(masses))]

    df = pd.DataFrame(list(zip(distances, masses, smoothed_masses, in_main_peak)),
                      columns=['distance', 'mass', 'smoothed_mass', 'in_main_peak'])

    g = (ggplot(df) +
         geom_segment(aes(x='distance', xend='distance', y=0, yend='mass', colour='in_main_peak'),
                      size=1) +
         scale_color_manual({True: '#7570b3', False: '#eeaaaa'}, guide=False) +
         geom_line(aes(x='distance', y='smoothed_mass'), size=0.5) +
         geom_vline(xintercept=mean, colour='#d95f02', linetype='dotted', size=0.5) +
         geom_vline(xintercept=peak_mean, colour='#d95f02', linetype='dashed', size=0.5) +
         theme_bw() +
         labs(title=title, x='distance', y='probability mass'))

    if sqrt_x:
        g += scale_x_sqrt(limits=(0, x_max))
    else:
        g += scale_x_continuous(limits=(0, x_max))

    if sqrt_y:
        g += scale_y_sqrt(expand=(0, 0), limits=(0, y_max))
    else:
        g += scale_y_continuous(expand=(0, 0), limits=(0, y_max))

    g.draw(show=True)


def load_distance_distribution(alignment_results, assembly_1, assembly_2):
    with open(alignment_results, 'rt') as results:
        for line in results:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            if parts[0] == assembly_1 and parts[1] == assembly_2:
                piece_size = int(parts[2])
                aligned_frac = float(parts[3])
                masses = [float(p) for p in parts[4:]]
                return piece_size, aligned_frac, masses
    sys.exit(f'\nError: could not find {assembly_1} and {assembly_2} in {alignment_results}')
