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
    scale_x_continuous, scale_x_sqrt, scale_y_continuous, scale_y_sqrt
import sys

from .distance import get_distance, get_peak_distance
from .log import log


def view(args):
    piece_size, aligned_frac, masses = \
        load_distance_distribution(args.alignment_results, args.assembly_1, args.assembly_2)

    title = f'{args.assembly_1} vs {args.assembly_2} ({piece_size} bp windows)'
    mean = get_distance(masses, piece_size, 'mean')
    peak, smoothed_masses = get_peak_distance(masses)
    peak /= piece_size
    log(f'Peak distance: {peak * piece_size}')  # TEMP

    x_max = len(masses) / piece_size
    y_max = 1.05 * max(max(masses), max(smoothed_masses))

    distances = [i / piece_size for i in range(len(masses))]
    df = pd.DataFrame(list(zip(distances, masses, smoothed_masses)),
                      columns=['distance', 'mass', 'smoothed_mass'])

    g = (ggplot(df) +
         geom_segment(aes(x='distance', xend='distance', y=0, yend='mass'),
                      size=1, colour='#7570b3') +
         geom_line(aes(x='distance', y='smoothed_mass'), size=0.5) +
         geom_vline(xintercept=mean, colour='#d95f02', linetype='dotted', size=0.5) +
         geom_vline(xintercept=peak, colour='#d95f02', linetype='dashed', size=0.5) +
         theme_bw() +
         labs(title=title))

    if args.sqrt_x:
        g += scale_x_sqrt(limits=(0, x_max))
    else:
        g += scale_x_continuous(limits=(0, x_max))

    if args.sqrt_y:
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
