"""
This module contains code for the 'verticall mask' subcommand.

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

import svgwrite
import sys

from .log import log, section_header, explanation, warning
from .misc import iterate_fasta, list_differences
from .tsv import get_column_index, check_header_for_assembly_a_regions, get_start_end


def mask(args):
    welcome_message(args)
    data, ref_name, ref_length, sample_names = load_regions(args.in_tsv, args.reference, args.multi)
    sequences, sample_names = load_pseudo_alignment(args.in_alignment, ref_name, sample_names)
    masked_sequences = mask_sequences(data, sequences, ref_name, ref_length, sample_names,
                                      args.h_char, args.u_char, args.image, args.vertical_colour,
                                      args.horizontal_colour, args.unaligned_colour)
    masked_sequences = finalise(masked_sequences, args.exclude_invariant)
    save_to_file(masked_sequences, args.out_alignment)
    finished_message()


def welcome_message(args):
    section_header('Starting Verticall mask')
    explanation('Verticall mask uses the results from Verticall pairwise to mask horizontal '
                'regions from a whole-genome pseudo-alignment.')
    log(f'Input pairwise tsv:      {args.in_tsv}')
    log(f'Input pseudo-alignment:  {args.in_alignment}')
    log(f'Output masked alignment: {args.out_alignment}')
    log()


def finished_message():
    section_header('Finished!')
    explanation('You can now use the masked pseudo-alignment to build a phylogeny using tools '
                'such as RAxML.')


def load_regions(filename, ref_name, multi):
    section_header('Loading input data')
    log(f'{filename}:')
    ref_name = check_tsv_file(filename, ref_name)
    data, distances, sample_names = {}, {}, set()
    v_col, h_col, u_col = None, None, None
    loaded_count = 0
    excluded_samples = get_multi_result_samples(filename, ref_name) if multi == 'exclude' else set()
    with open(filename, 'rt') as pairwise_file:
        for i, line in enumerate(pairwise_file):
            parts = line.strip('\n').split('\t')
            if i == 0:
                v_col = get_column_index(parts, 'assembly_a_vertical_regions', filename)
                h_col = get_column_index(parts, 'assembly_a_horizontal_regions', filename)
                u_col = get_column_index(parts, 'assembly_a_unaligned_regions', filename)
                d_col = get_column_index(parts, 'mean_vertical_distance', filename)
                continue
            if parts[0] != ref_name:
                continue
            assembly_name = parts[1]
            if assembly_name in excluded_samples:
                continue
            distance = float(parts[d_col])
            regions = load_regions_one_assembly(parts, v_col, h_col, u_col)
            if assembly_name not in data:  # first time we've seen this assembly
                distances[assembly_name] = distance
                data[assembly_name] = regions
                loaded_count += 1
                sample_names.add(assembly_name)
            else:  # seen this assembly already
                if multi == 'first':
                    pass
                elif multi == 'low':
                    if distance < distances[assembly_name]:
                        distances[assembly_name] = distance
                        data[assembly_name] = regions
                elif multi == 'high':
                    if distance > distances[assembly_name]:
                        distances[assembly_name] = distance
                        data[assembly_name] = regions
                else:
                    assert False
    if len(data) == 0:
        sys.exit(f'Error: no reference-to-assembly pairwise comparisons found in {filename} - is '
                 f'the provided reference name correct?')
    log(f'  {loaded_count} reference-to-assembly pairwise comparisons')
    ref_length = get_ref_length(data)
    log()
    return data, ref_name, ref_length, sorted(sample_names)


def get_multi_result_samples(filename, ref_name):
    samples, multi_result_samples = set(), set()
    with open(filename, 'rt') as f:
        for i, line in enumerate(f):
            if i == 0:  # header line
                continue
            parts = line.strip('\n').split('\t')
            assembly_a, assembly_b = parts[0], parts[1]
            if assembly_a != ref_name:
                continue
            if assembly_b in samples:
                multi_result_samples.add(assembly_b)
            else:
                samples.add(assembly_b)
    if multi_result_samples:
        multi_result_samples_str = ', '.join(sorted(multi_result_samples))
        warning(f'The following samples will be excluded due to secondary results: '
                f'{multi_result_samples_str}')
    return multi_result_samples


def check_tsv_file(filename, ref_name):
    a_sample_names = set()
    with open(filename, 'rt') as f:
        for i, line in enumerate(f):
            parts = line.strip('\n').split('\t')
            if i == 0:  # header line
                check_header_for_assembly_a_regions(parts, filename)
            else:
                a_sample_names.add(parts[0])
    if ref_name is None:  # user didn't specify a reference name
        if len(a_sample_names) == 1:
            ref_name = list(a_sample_names)[0]
            log(f'  Automatically determined reference name: {ref_name}')
        else:
            sys.exit('Error: could not automatically determine the reference name, please '
                     'specify one using the --reference option.')
    return ref_name


def load_regions_one_assembly(parts, v_col, h_col, u_col):
    contig_names = set()

    vertical_regions = parts[v_col].split(',') if parts[v_col] else []
    contig_names |= {r.split(':')[0] for r in vertical_regions}
    vertical_regions = [get_start_end(r) for r in vertical_regions]

    horizontal_regions = parts[h_col].split(',') if parts[h_col] else []
    contig_names |= {r.split(':')[0] for r in horizontal_regions}
    horizontal_regions = [get_start_end(r) for r in horizontal_regions]

    unaligned_regions = parts[u_col].split(',') if parts[u_col] else []
    contig_names |= {r.split(':')[0] for r in unaligned_regions}
    unaligned_regions = [get_start_end(r) for r in unaligned_regions]

    if len(contig_names) > 1:
        contig_names_str = ', '.join(sorted(contig_names))
        sys.exit(f'Error: reference genome has more than one contig name ({contig_names_str})')

    # Double check that the data makes sense - the entire reference sequence should be covered once.
    all_regions = sorted(vertical_regions + horizontal_regions + unaligned_regions)
    for i, r in enumerate(all_regions):
        start, end = r
        if i > 0:
            prev_start, prev_end = all_regions[i-1]
            assert start == prev_end

    return vertical_regions, horizontal_regions, unaligned_regions


def get_ref_length(data):
    ref_lengths = set()
    for vertical_regions, horizontal_regions, unaligned_regions in data.values():
        ref_length = max(r[1] for r in vertical_regions + horizontal_regions + unaligned_regions)
        ref_lengths.add(ref_length)
    if len(ref_lengths) > 1:
        sys.exit('Error: multiple inconsistent reference sequence lengths')
    ref_length = list(ref_lengths)[0]
    log(f'  Reference length: {ref_length:,} bp')
    return ref_length


def load_pseudo_alignment(filename, ref_name, tsv_sample_names):
    alignment = {}
    log(f'{filename}:')
    sequence_names, sequence_lengths = set(), set()
    for name, seq in iterate_fasta(filename, preserve_case=True):
        sequence_names.add(name)
        sequence_lengths.add(len(seq))
        alignment[name] = seq

    if len(alignment) == 0:
        sys.exit(f'Error: no sequences could be loaded from {filename}')
    log(f'  {len(sequence_names)} sequences loaded')

    if len(sequence_lengths) != 1:
        sys.exit(f'Error: all sequences in {filename} must be the same length')
    sequence_length = list(sequence_lengths)[0]
    log(f'  Pseudo-alignment length: {sequence_length:,} bp')
    log()

    if ref_name not in sequence_names:
        sys.exit(f'Error: could not find reference sequence ({ref_name}) in {filename}')
    sequence_names.discard(ref_name)

    sequence_names = sorted(sequence_names)
    in_both, in_tsv_not_alignment, in_alignment_not_tsv = list_differences(tsv_sample_names,
                                                                           sequence_names)
    if len(in_both) == 0:
        sys.exit('Error: tsv and pseudo-alignment files have no sample names in common')
    if len(in_tsv_not_alignment) == 0 and len(in_alignment_not_tsv) == 0:
        log(f'All {len(in_both)} sample names match in tsv and pseudo-alignment files')
    else:
        log(f'{len(in_both)} sample names are common to both tsv and pseudo-alignment files')
    log()

    return alignment, in_both


def mask_sequences(data, sequences, ref_name, ref_length, sample_names, h_char, u_char,
                   image_filename, v_colour, h_colour, u_colour):
    section_header('Masking sequences')
    ref_seq = sequences[ref_name]
    ref_pos_to_align_pos = get_alignment_positions(ref_seq, ref_length)
    longest_sample_name_len = max(len(s) for s in sample_names)
    masked_sequences = {ref_name: ref_seq}

    if image_filename is not None:
        image = svgwrite.Drawing(image_filename, profile='full')
    else:
        image = None

    y_pos = 12
    for sample_name in sample_names:
        log(f'{sample_name.rjust(longest_sample_name_len)}:', end=' ')
        masked_sequences[sample_name] = \
            mask_one_sequence(data, sequences, sample_name, h_char, u_char, ref_pos_to_align_pos,
                              ref_length, image, v_colour, h_colour, u_colour, y_pos)
        y_pos += 12

    if image_filename is not None:
        image.save()
    log()
    return masked_sequences


def mask_one_sequence(data, sequences, sample_name, h_char, u_char, ref_pos_to_align_pos,
                      ref_length, image, v_colour, h_colour, u_colour, y_pos):
    _, horizontal_regions, unaligned_regions = data[sample_name]
    sample_seq = [b for b in sequences[sample_name]]
    unmasked, h_masked, u_masked = ref_length, 0, 0
    if image is not None:
        image.add(image.text(sample_name, insert=(97, y_pos+4), style='text-anchor:end',
                             font_size='12px'))
        image.add(image.line((100, y_pos), (500, y_pos), stroke=v_colour, stroke_width=9))
    if h_char is not None:
        for start, end in horizontal_regions:
            unmasked -= (end - start)
            h_masked += (end - start)
            start, end = ref_pos_to_align_pos[start], ref_pos_to_align_pos[end]
            for i in range(start, end):
                sample_seq[i] = h_char
            if image is not None:
                image.add(image.line((100 + 400 * start / ref_length, y_pos),
                                     (100 + 400 * end / ref_length, y_pos),
                                     stroke=h_colour, stroke_width=9))
    if u_char is not None:
        for start, end in unaligned_regions:
            unmasked -= (end - start)
            u_masked += (end - start)
            start, end = ref_pos_to_align_pos[start], ref_pos_to_align_pos[end]
            for i in range(start, end):
                sample_seq[i] = u_char
            if image is not None:
                image.add(image.line((100 + 400 * start / ref_length, y_pos),
                                     (100 + 400 * end / ref_length, y_pos),
                                     stroke=u_colour, stroke_width=9))
    log_message = f'{100.0 * unmasked/ref_length:6.2f}% unmasked'
    if h_char is not None:
        log_message += f', {100.0 * h_masked/ref_length:5.2f}% "{h_char}"'
    if u_char is not None:
        log_message += f', {100.0 * u_masked/ref_length:5.2f}% "{u_char}"'
    log(log_message)
    return ''.join(sample_seq)


def get_alignment_positions(aligned_ref_seq, ref_length):
    """
    Returns a dictionary that translates reference positions to alignment positions. If the
    alignment contains no insertions in the reference sequence, these two sets of positions will be
    the same, but if there are insertions in the reference sequence, then the alignment positions
    can be bigger than the reference positions.
    """
    positions = {}
    ref_pos = 0
    for alignment_pos, base in enumerate(aligned_ref_seq):
        positions[ref_pos] = alignment_pos
        if base != '-':
            ref_pos += 1
    positions[ref_pos] = len(aligned_ref_seq)
    if ref_pos != ref_length:
        sys.exit('Error: length of reference sequence in alignment does not match length of '
                 'reference sequence in TSV file - have regions been masked with dashes?')
    return positions


def finalise(masked_sequences, exclude_invariant):
    section_header('Finalising sequences')
    masked_sequences = drop_gap_positions(masked_sequences)
    if exclude_invariant:
        masked_sequences = drop_invariant_positions(masked_sequences)
    log()
    return masked_sequences


def save_to_file(masked_sequences, filename):
    log(f'Saving masked pseudo-alignment to {filename}')
    with open(filename, 'wt') as f:
        for name, seq in masked_sequences.items():
            if len(seq) == 0:
                warning(f'excluded {name} due to empty sequence')
            else:
                f.write(f'>{name}\n{seq}\n')
    log()


def drop_gap_positions(sequences):
    """
    Returns an alignment where any columns that consist entirely of gaps are removed
    """
    alignment_length = get_alignment_length(sequences)
    positions_to_remove = set()
    for i in range(alignment_length):
        for seq in sequences.values():
            base = seq[i].upper()
            if base != '-':
                break
        else:  # no break
            positions_to_remove.add(i)
    log(f'{len(positions_to_remove):,} all-gap positions '
        f'({100.0 * len(positions_to_remove)/alignment_length:.3}%) removed from pseudo-alignment')
    return drop_positions(sequences, positions_to_remove)


def drop_invariant_positions(sequences):
    """
    Returns an alignment where any columns that lack variation are removed.
    """
    alignment_length = get_alignment_length(sequences)
    positions_to_remove = set()
    for i in range(alignment_length):
        bases_at_pos = {seq[i] for seq in sequences.values()}
        if count_real_bases(bases_at_pos) < 2:
            positions_to_remove.add(i)
    log(f'{len(positions_to_remove):,} invariant positions '
        f'({100.0 * len(positions_to_remove)/alignment_length:.3}%) removed from pseudo-alignment')
    # TODO: display what the invariant positions were (counts for A, C, G, T and ambiguous)
    return drop_positions(sequences, positions_to_remove)


def get_alignment_length(sequences):
    alignment_lengths = {len(seq) for seq in sequences.values()}
    assert len(alignment_lengths) == 1
    return list(alignment_lengths)[0]


def drop_positions(sequences, positions_to_remove):
    if len(positions_to_remove) == 0:
        return sequences
    new_sequences = {}
    for name, seq in sequences.items():
        new_seq = ''.join(b for i, b in enumerate(seq) if i not in positions_to_remove)
        new_sequences[name] = new_seq
    return new_sequences


def count_real_bases(base_set):
    count = 0
    if 'A' in base_set or 'a' in base_set:
        count += 1
    if 'C' in base_set or 'c' in base_set:
        count += 1
    if 'G' in base_set or 'g' in base_set:
        count += 1
    if 'T' in base_set or 't' in base_set:
        count += 1
    return count
