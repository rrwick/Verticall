#!/usr/bin/env python3
"""
This module contains code for the 'verticall pairwise' subcommand.

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

import collections
import gzip
import re

from .log import log
from .misc import iterate_fasta, contains_ambiguous_bases


def repair(args):
    log()
    log(f'Repairing {args.in_file}')

    new_names, new_infos, new_seqs = [], [], []
    before_contigs, before_bases = 0, 0
    for name, info, seq in iterate_fasta(args.in_file, include_info=True):
        before_contigs += 1
        before_bases += len(seq)
        for seq_part in split_seq_on_ambiguous(seq):
            new_names.append(name)
            new_infos.append(info)
            new_seqs.append(seq_part)

    contig = 'contig' if before_contigs == 1 else 'contigs'
    log(f'  before repair: {before_contigs} {contig}, {before_bases} bp')
    contig = 'contig' if len(new_names) == 1 else 'contigs'
    log(f'  after repair:  {len(new_names)} {contig}, {sum(len(s) for s in new_seqs)} bp')

    new_names = make_names_unique(new_names)
    if args.in_file == args.out_file:
        log(f'  overwriting original assembly file')
    else:
        log(f'  writing new assembly file: {args.out_file}')
    save_seqs_to_file(new_names, new_infos, new_seqs, args.out_file)


def save_seqs_to_file(names, infos, seqs, filename):
    assert len(names) == len(infos) == len(seqs)
    if str(filename).endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open
    with open_func(filename, 'wt') as f:
        for name, info, seq in zip(names, infos, seqs):
            if info:
                header = f'>{name} {info}'
            else:
                header = f'>{name}'
            f.write(f'{header}\n{seq}\n')


def split_seq_on_ambiguous(seq):
    seq_parts = re.split('[^ACGT]', seq.upper())
    return [p for p in seq_parts if p]


def make_names_unique(names):
    """
    Given a list of contig names, this function returns a new list of the same length where all
    names are unique. It does this by appending underscore+number for duplicates. It can also call
    itself recursively in case adding these suffixes creates new duplicates.
    """
    duplicates = {name for name, count in collections.Counter(names).items() if count > 1}
    times_used = collections.defaultdict(int)
    unique_names = []
    any_changes = False
    for name in names:
        if name in duplicates:
            times_used[name] += 1
            unique_names.append(name + '_' + str(times_used[name]))
            any_changes = True
        else:
            unique_names.append(name)
    if any_changes:
        return make_names_unique(unique_names)
    else:
        return unique_names
