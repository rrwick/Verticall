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

import numpy as np
import pandas as pd
from plotnine import ggplot, aes, geom_segment, geom_vline, geom_line, labs, theme_bw, \
    scale_y_continuous, scale_y_sqrt
from scipy.stats import gamma
import sys

from .distance import get_distance, get_tightest_half
from .gamma import fit_gamma_to_distribution


def view(args):
    piece_size, aligned_frac, masses = \
        load_distance_distribution(args.alignment_results, args.assembly_1, args.assembly_2)

    distances = [i / piece_size for i in range(len(masses))]
    df = pd.DataFrame(list(zip(distances, masses)),  columns=['distance', 'mass'])
    title = f'{args.assembly_1} vs {args.assembly_2}'
    mean = get_distance(masses, piece_size, 'mean')
    median = get_distance(masses, piece_size, 'median')
    median_int = get_distance(masses, piece_size, 'median_int')

    low, high = get_tightest_half(masses)
    low /= piece_size
    high /= piece_size

    shape, scale, vertical_scale = fit_gamma_to_distribution(masses)
    gamma_x = np.arange(len(masses))
    gamma_y = gamma.pdf(gamma_x, shape, scale=scale)
    gamma_x = gamma_x / piece_size
    gamma_y = gamma_y * vertical_scale
    gamma_df = pd.DataFrame({'x': gamma_x, 'y': gamma_y})

    g = (ggplot(df, aes('distance', 'mass')) +
         geom_segment(aes(x='distance', xend='distance', y=0, yend='mass'),
                      colour='#880000', size=1) +
         geom_line(data=gamma_df, mapping=aes(x='x', y='y')) +
         geom_vline(xintercept=mean, colour='#008888', linetype='dotted') +
         geom_vline(xintercept=median, colour='#0000bb', linetype='dotted') +
         geom_vline(xintercept=median_int, colour='#00bb00', linetype='dotted') +
         geom_vline(xintercept=low, colour='#aaaaaa', linetype='dashed') +
         geom_vline(xintercept=high, colour='#aaaaaa', linetype='dashed') +
         theme_bw() +
         labs(title=title))

    y_max = 1.05 * max(masses)
    if args.sqrt_y:
        g += scale_y_sqrt(expand=(0, 0), limits=(0, y_max))
    else:
        g += scale_y_continuous(expand=(0, 0), limits=(0, y_max))

    print(g)


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
