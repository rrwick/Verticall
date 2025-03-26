"""
This module contains code for the 'verticall view' subcommand.

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

import pandas as pd
import warnings
import sys

from .distance import get_distance
from .log import log, section_header, explanation
from .pairwise import find_assemblies, check_assemblies, build_indices, process_one_pair, \
    prepare_log_text


warnings.filterwarnings('ignore')


def view(args):
    welcome_message()
    assemblies = find_assemblies(args.in_dir)
    name_a, name_b, filename_a, filename_b = get_sample_names_and_filenames(args, assemblies)
    assemblies = [(name_a, filename_a), (name_b, filename_b)]
    check_assemblies(assemblies)
    build_indices(args, assemblies)
    section_header('Processing assembly pair')
    all_args = args, name_a, name_b, filename_a, filename_b
    alignments, window_size, masses, smoothed_masses, thresholds, vertical_masses, \
        horizontal_masses, painted_a, log_text = process_one_pair(all_args, view=True,
                                                                  view_num=args.result)
    log()
    log('\n'.join(prepare_log_text(log_text, True)), end='\n\n')
    finished_message()
    show_plots(name_a, name_b, alignments, window_size, masses, smoothed_masses, thresholds,
               vertical_masses, horizontal_masses, painted_a, args.sqrt_distance, args.sqrt_mass,
               args.vertical_colour, args.horizontal_colour, args.ambiguous_colour)


def welcome_message():
    section_header('Starting Verticall view')
    explanation('Verticall view performs the same analysis as Verticall pairwise, but only for a '
                'single assembly pair. The results of the analysis are then visualised in '
                'interactive plots.')


def finished_message():
    section_header('Finished!')
    explanation('The plot will now be displayed. To complete the process, close the interactive '
                'plot windows.')


def get_sample_names_and_filenames(args, assemblies):
    name_a, name_b = args.names.split(',')
    try:
        filename_a = [filename for name, filename in assemblies if name == name_a][0]
    except IndexError:
        sys.exit(f'Error: could not find assembly named {name_a}')
    try:
        filename_b = [filename for name, filename in assemblies if name == name_b][0]
    except IndexError:
        sys.exit(f'Error: could not find assembly named {name_b}')
    return name_a, name_b, filename_a, filename_b


def show_plots(sample_name_a, sample_name_b, alignments, window_size, masses, smoothed_masses,
               thresholds, vertical_masses, horizontal_masses, painted_a, sqrt_distance,
               sqrt_mass, vertical_colour, horizontal_colour, ambiguous_colour):
    import matplotlib.pyplot as plt
    fig_1 = distribution_plot_1(sample_name_a, sample_name_b, window_size, masses, smoothed_masses,
                                thresholds, sqrt_distance, sqrt_mass, vertical_colour,
                                horizontal_colour, ambiguous_colour)
    fig_2 = distribution_plot_2(sample_name_a, sample_name_b, window_size, vertical_masses,
                                horizontal_masses, sqrt_distance, sqrt_mass, vertical_colour,
                                horizontal_colour)
    fig_3 = alignment_plot(sample_name_a, sample_name_b, alignments, window_size, sqrt_distance,
                           vertical_colour, horizontal_colour, ambiguous_colour)
    fig_4 = contig_plot(sample_name_a, painted_a, window_size, sqrt_distance, vertical_colour,
                        horizontal_colour)

    plt.show()


def distribution_plot_1(sample_name_a, sample_name_b, window_size, masses, smoothed_masses,
                        thresholds, sqrt_distance, sqrt_mass, vertical_colour, horizontal_colour,
                        ambiguous_colour):
    from plotnine import ggplot, aes, geom_segment, geom_line, geom_vline, labs, theme_bw, \
        scale_x_continuous, scale_x_sqrt, scale_y_continuous, scale_y_sqrt, scale_color_manual, \
        theme
    title = f'{sample_name_a} vs {sample_name_b} full distribution with thresholds'
    mean = get_distance(masses, window_size, 'mean')
    median = get_distance(masses, window_size, 'median')
    x_max = len(masses) / window_size
    y_max = 1.05 * max(max(masses), max(smoothed_masses))
    distances = [i / window_size for i in range(len(masses))]
    grouping = group_using_thresholds(masses, thresholds)

    df = pd.DataFrame(list(zip(distances, masses, smoothed_masses, grouping)),
                      columns=['distance', 'mass', 'smoothed_mass', 'grouping'])

    g = (ggplot(df) +
         geom_segment(aes(x='distance', xend='distance', y=0, yend='mass', colour='grouping'),
                      size=1) +
         scale_color_manual({'very_low': horizontal_colour, 'low': ambiguous_colour,
                             'very_high': horizontal_colour, 'high': ambiguous_colour,
                             'central': vertical_colour}, guide=None) +
         geom_line(aes(x='distance', y='smoothed_mass'), size=0.5) +
         geom_vline(xintercept=mean, colour='#000000', linetype='dotted', size=0.5) +
         geom_vline(xintercept=median, colour='#000000', linetype='dashed', size=0.5) +
         theme_bw() + theme(figure_size=(10, 5)) +
         labs(title=title, x='distance', y='mass'))

    if sqrt_distance:
        g += scale_x_sqrt(limits=(0, x_max))
    else:
        g += scale_x_continuous(limits=(0, x_max))
    if sqrt_mass:
        g += scale_y_sqrt(expand=(0, 0), limits=(0, y_max))
    else:
        g += scale_y_continuous(expand=(0, 0), limits=(0, y_max))

    return g.draw(show=True)


def distribution_plot_2(sample_name_a, sample_name_b, window_size, vertical_masses,
                        horizontal_masses, sqrt_distance, sqrt_mass, vertical_colour,
                        horizontal_colour):
    from plotnine import ggplot, aes, geom_segment, geom_vline, labs, theme_bw, \
        scale_x_continuous, scale_x_sqrt, scale_y_continuous, scale_y_sqrt, theme
    title = f'{sample_name_a} vs {sample_name_b} vertical vs horizontal distribution'
    mean = get_distance(vertical_masses, window_size, 'mean')
    median = get_distance(vertical_masses, window_size, 'median')
    max_distance = max(len(vertical_masses), len(horizontal_masses))
    x_max = max_distance / window_size
    y_max = 1.05 * max(max(vertical_masses), max(horizontal_masses))

    distances = [i / window_size for i in range(max_distance)]

    total_masses = [vertical_masses[i] + horizontal_masses[i] for i in range(max_distance)]

    df = pd.DataFrame(list(zip(distances, vertical_masses, horizontal_masses, total_masses)),
                      columns=['distance', 'vertical_mass', 'horizontal_mass', 'total_mass'])

    g = (ggplot(df) +
         geom_segment(aes(x='distance', xend='distance', y=0, yend='vertical_mass'),
                      size=1, colour=vertical_colour) +
         geom_segment(aes(x='distance', xend='distance', y='vertical_mass', yend='total_mass'),
                      size=1, colour=horizontal_colour) +
         geom_vline(xintercept=mean, colour='#000000', linetype='dotted', size=0.5) +
         geom_vline(xintercept=median, colour='#000000', linetype='dashed', size=0.5) +
         theme_bw() + theme(figure_size=(10, 5)) +
         labs(title=title, x='distance', y='mass'))

    if sqrt_distance:
        g += scale_x_sqrt(limits=(0, x_max))
    else:
        g += scale_x_continuous(limits=(0, x_max))
    if sqrt_mass:
        g += scale_y_sqrt(expand=(0, 0), limits=(0, y_max))
    else:
        g += scale_y_continuous(expand=(0, 0), limits=(0, y_max))

    return g.draw(show=True)


def group_using_thresholds(masses, thresholds):
    """
    Assigns a group (very_low, low, central, high, very_high) to each mass, returned as a list.
    """
    grouping = []
    very_low, low = thresholds['very_low'], thresholds['low']
    very_high, high = thresholds['very_high'], thresholds['high']
    for i in range(len(masses)):
        if very_low is not None and i < very_low:
            grouping.append('very_low')
        elif low is not None and i < low:
            grouping.append('low')
        elif very_high is not None and i > very_high:
            grouping.append('very_high')
        elif high is not None and i > high:
            grouping.append('high')
        else:
            grouping.append('central')
    return grouping


def alignment_plot(sample_name_a, sample_name_b, alignments, window_size, sqrt_distance,
                   vertical_colour, horizontal_colour, ambiguous_colour, include_ambiguous=False):
    from plotnine import ggplot, aes, geom_line, geom_vline, labs, theme_bw, scale_x_continuous, \
        scale_y_continuous, scale_y_sqrt, element_blank, theme, annotate
    title = f'{sample_name_a} vs {sample_name_b} painted alignments'

    boundaries = [0]
    x_max = 0
    max_differences = 1
    for a in alignments:
        x_max += len(a.simplified_cigar)
        boundaries.append(x_max)
        max_differences = max(max_differences, a.get_max_differences())
    y_max = 1.05 * (max_differences / window_size)

    g = (ggplot() +
         theme_bw() + theme(figure_size=(15, 5)) +
         theme(panel_grid_major_x=element_blank(), panel_grid_minor_x=element_blank()) +
         scale_x_continuous(expand=(0, 0), limits=(0, x_max)) +
         labs(title=title, x='alignment position', y='distance'))

    if sqrt_distance:
        g += scale_y_sqrt(expand=(0, 0), limits=(0, y_max))
    else:
        g += scale_y_continuous(expand=(0, 0), limits=(0, y_max))

    for b in boundaries:
        g += geom_vline(xintercept=b, colour='#555555', size=0.5)

    offset = 0
    for a in alignments:
        for start, end in a.get_vertical_blocks(include_ambiguous):
            g += annotate('rect', xmin=start+offset, xmax=end+offset, ymin=0.0, ymax=y_max,
                          fill=vertical_colour, alpha=0.35)
        for start, end in a.get_horizontal_blocks(include_ambiguous):
            g += annotate('rect', xmin=start+offset, xmax=end+offset, ymin=0.0, ymax=y_max,
                          fill=horizontal_colour, alpha=0.35)
        for start, end in a.get_ambiguous_blocks(include_ambiguous):
            g += annotate('rect', xmin=start+offset, xmax=end+offset, ymin=0.0, ymax=y_max,
                          fill=ambiguous_colour, alpha=0.35)
        positions = [offset + ((w[0] + w[1]) / 2.0) for w in a.windows_no_overlap]
        distances = [d / window_size for d in a.window_differences]
        df = pd.DataFrame(list(zip(positions, distances)), columns=['pos', 'dist'])
        g += geom_line(data=df, mapping=aes(x='pos', y='dist'), size=0.5)
        offset += len(a.simplified_cigar)

    return g.draw(show=True)


def contig_plot(sample_name, painted, window_size, sqrt_distance, vertical_colour,
                horizontal_colour):
    from plotnine import ggplot, aes, geom_line, geom_vline, labs, theme_bw, scale_x_continuous, \
        scale_y_continuous, scale_y_sqrt, element_blank, theme, annotate
    title = f'{sample_name} painted contigs'

    boundaries = [0]
    x_max = 0
    for name, contig in painted.contigs.items():
        x_max += contig.length
        boundaries.append(x_max)
    y_max = 1.05 * (max(1, painted.get_max_differences()) / window_size)

    g = (ggplot() +
         theme_bw() + theme(figure_size=(15, 5)) +
         theme(panel_grid_major_x=element_blank(), panel_grid_minor_x=element_blank()) +
         scale_x_continuous(expand=(0, 0), limits=(0, x_max)) +
         labs(title=title, x='contig position', y='distance'))

    if sqrt_distance:
        g += scale_y_sqrt(expand=(0, 0), limits=(0, y_max))
    else:
        g += scale_y_continuous(expand=(0, 0), limits=(0, y_max))

    for b in boundaries:
        g += geom_vline(xintercept=b, colour='#555555', size=0.5)

    offset = 0
    for name, contig in painted.contigs.items():
        for start, end in contig.get_vertical_blocks():
            g += annotate('rect', xmin=start+offset, xmax=end+offset, ymin=0.0, ymax=y_max,
                          fill=vertical_colour, alpha=0.35)
        for start, end in contig.get_horizontal_blocks():
            g += annotate('rect', xmin=start+offset, xmax=end+offset, ymin=0.0, ymax=y_max,
                          fill=horizontal_colour, alpha=0.35)
        for points in contig.alignment_points:
            positions = [offset + p[0] for p in points]
            distances = [p[1] / window_size for p in points]
            df = pd.DataFrame(list(zip(positions, distances)), columns=['pos', 'dist'])
            g += geom_line(data=df, mapping=aes(x='pos', y='dist'), size=0.5)
        offset += contig.length

    return g.draw(show=True)


def check_hex_colour(colour):
    if len(colour) != 7:
        return False
    if colour[0] != '#':
        return False
    colour = colour.lower()
    good_chars = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'}
    for i in range(1, 7):
        if colour[i] not in good_chars:
            return False
    return True
