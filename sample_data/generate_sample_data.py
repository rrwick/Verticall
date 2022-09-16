#!/usr/bin/env python3
"""
This is a standalone script to generate synthetic sample data for Verticall. It starts with a
random sequence (used as the reference) and 'evolves' it into a simple tree of five isolates with a
bit of recombination. It produces an 'assemblies' directory containing five FASTA files, a
'reference.fasta' file and an 'alignment.fasta' file.

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

import pathlib
import random

SEQ_LENGTH = 2000000
RECOMBINATION = True


def main():
    random.seed(0)

    #                             ↓
    #                   ┌─────────T─────────A
    #           ┌───────Y
    #           │       └───────────────U───B
    #           │                       ↑
    # R─────────Z                   ┌───V───C
    #           │           ┌───────W
    #           └───────────X       └───────D
    #                       │
    #                       └───────────────E
    r = get_random_seq(SEQ_LENGTH)
    z = mutate_seq(r, 0.005)
    y = mutate_seq(z, 0.004)
    t = mutate_seq(y, 0.005)
    if RECOMBINATION:
        start, end = int(2*SEQ_LENGTH/8), int(3*SEQ_LENGTH/8)
        t = t[0:start] + mutate_seq(r[start:end], 0.05) + t[end:]
    u = mutate_seq(y, 0.008)
    x = mutate_seq(z, 0.006)
    w = mutate_seq(x, 0.004)
    v = mutate_seq(w, 0.002)
    if RECOMBINATION:
        start, end = int(12*SEQ_LENGTH/16), int(13*SEQ_LENGTH/16)
        u = u[0:start] + v[start:end] + u[end:]
    a = mutate_seq(t, 0.005)
    b = mutate_seq(u, 0.002)
    c = mutate_seq(v, 0.002)
    d = mutate_seq(w, 0.004)
    e = mutate_seq(x, 0.008)

    save_to_fasta('alignment.fasta',
                  [('reference', r), ('a', a), ('b', b), ('c', c), ('d', d), ('e', e)])
    save_to_fasta('reference.fasta', [('reference', r)])
    pathlib.Path('assemblies').mkdir(exist_ok=True)
    save_to_fasta('assemblies/a.fasta', [('a', rotate_seq(a))])
    save_to_fasta('assemblies/b.fasta', [('b', rotate_seq(b))])
    save_to_fasta('assemblies/c.fasta', [('c', rotate_seq(c))])
    save_to_fasta('assemblies/d.fasta', [('d', rotate_seq(d))])
    save_to_fasta('assemblies/e.fasta', [('e', rotate_seq(e))])


def get_random_base():
    return {0: 'A', 1: 'C', 2: 'G', 3: 'T'}[random.randint(0, 3)]


def get_random_different_base(b):
    random_base = get_random_base()
    while b == random_base:
        random_base = get_random_base()
    return random_base


def get_random_seq(seq_len):
    return ''.join(get_random_base() for _ in range(seq_len))


def mutate_seq(seq, distance):
    new_seq = []
    for base in seq:
        if random.random() < distance:
            new_seq.append(get_random_different_base(base))
        else:
            new_seq.append(base)
    return ''.join(new_seq)


def save_to_fasta(filename, seqs):
    with open(filename, 'wt') as f:
        for name, seq in seqs:
            f.write(f'>{name}\n{seq}\n')


def reverse_complement(seq):
    rev_bases = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join([rev_bases[x] for x in seq][::-1])


def rotate_seq(seq):
    if random.choice([0, 1]) == 0:
        seq = reverse_complement(seq)
    pos = random.randint(0, len(seq) - 1)
    return seq[pos:] + seq[:pos]


if __name__ == '__main__':
    main()
