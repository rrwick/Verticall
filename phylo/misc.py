"""
Copyright 2022 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/XXXXXXXXX

This file is part of XXXXXXXXX. XXXXXXXXX is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. XXXXXXXXX is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with XXXXXXXXX.
If not, see <http://www.gnu.org/licenses/>.
"""

import gzip
import os
import subprocess
import sys

from .log import log, bold_yellow


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    unknown_file = open(str(filename), 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('\nError: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('\nError: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def check_python_version():
    if sys.version_info.major < 3 or sys.version_info.minor < 7:
        sys.exit('\nError: XXXXXXXXX requires Python 3.7 or later')


def get_ascii_art():
    ascii_art = (bold_yellow(r" _____   _             _        ") + '\n' +
                 bold_yellow(r"|  __ \ | |           | |       ") + '\n' +
                 bold_yellow(r"| |__) || |__   _   _ | |  ___  ") + '\n' +
                 bold_yellow(r"|  ___/ | '_ \ | | | || | / _ \ ") + '\n' +
                 bold_yellow(r"| |     | | | || |_| || || (_) |") + '\n' +
                 bold_yellow(r"|_|     |_| |_| \__, ||_| \___/ ") + '\n' +
                 bold_yellow(r"                 __/ |          ") + '\n' +
                 bold_yellow(r"                |___/         ") + '\n')
    return ascii_art


def get_sequence_file_type(filename):
    """
    Peeks into a file and returns either 'FASTA', 'FASTQ' or 'GFA' based on what the start of the
    file looks like.
    """
    if not os.path.isfile(filename):
        sys.exit(f'\nError: could not find {filename}')
    with get_open_func(filename)(filename, 'rt') as seq_file:
        try:
            first_char = seq_file.read(1)
        except UnicodeDecodeError:
            first_char = ''
    if first_char == '>':
        return 'FASTA'
    elif first_char == '@':
        return 'FASTQ'
    elif first_char == 'H' or first_char == 'S' or first_char == 'L':
        return 'GFA'
    else:
        return 'unknown'


def iterate_fasta(filename):
    """
    Takes a FASTA file as input and yields the contents as (name, seq) tuples.
    """
    with get_open_func(filename)(filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    yield name.split()[0], ''.join(sequence)
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line.upper())
        if name:
            yield name.split()[0], ''.join(sequence)


def iterate_gfa(filename):
    """
    Takes a GFA file as input and yields the contents as (name, seq) tuples. Only 'S' lines are
    used - everything else in the GFA file is ignored.
    """
    with get_open_func(filename)(filename, 'rt') as gfa_file:
        for line in gfa_file:
            if not line.startswith('S\t'):
                continue
            parts = line.strip().split('\t')
            yield parts[1], parts[2]


def get_generator(filename):
    """
    Returns either the iterate_fasta generator or the iterate_gfa generator, as appropriate for the
    file type.
    """
    file_type = get_sequence_file_type(filename)
    if file_type == 'FASTA':
        return iterate_fasta
    elif file_type == 'GFA':
        return iterate_gfa


def check_assembly_file_type(a):
    file_type = get_sequence_file_type(a)
    if file_type != 'FASTA' and file_type != 'GFA':
        sys.exit(f'\nError: {a} is not in FASTA or GFA format')
