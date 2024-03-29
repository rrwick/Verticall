"""
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

import gzip
import multiprocessing
import os
import sys

from .log import bold_yellow


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    https://stackoverflow.com/questions/13044562
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
        sys.exit('\nError: Verticall requires Python 3.7 or later')


def get_ascii_art():
    ascii_art = (bold_yellow(r" __      __          _    _               _  _ ") + '\n' +
                 bold_yellow(r" \ \    / /         | |  (_)             | || |") + '\n' +
                 bold_yellow(r"  \ \  / /___  _ __ | |_  _   ___   __ _ | || |") + '\n' +
                 bold_yellow(r"   \ \/ // _ \| '__|| __|| | / __| / _` || || |") + '\n' +
                 bold_yellow(r"    \  /|  __/| |   | |_ | || (__ | (_| || || |") + '\n' +
                 bold_yellow(r"     \/  \___||_|    \__||_| \___| \__,_||_||_|") + '\n')
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


def iterate_fasta(filename, include_info=False, preserve_case=False):
    """
    Takes a FASTA file as input and yields the contents as (name, seq) tuples. If include_info is
    set, it will yield (name, info, seq) tuples, where info is whatever follows the name.
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
                    if include_info:
                        name_parts = name.split(maxsplit=1)
                        contig_name = name_parts[0]
                        info = '' if len(name_parts) == 1 else name_parts[1]
                        yield contig_name, info, ''.join(sequence)
                    else:
                        yield name.split()[0], ''.join(sequence)
                    sequence = []
                name = line[1:]
            else:
                if preserve_case:
                    sequence.append(line)
                else:
                    sequence.append(line.upper())
        if name:
            if include_info:
                name_parts = name.split(maxsplit=1)
                contig_name = name_parts[0]
                info = '' if len(name_parts) == 1 else name_parts[1]
                yield contig_name, info, ''.join(sequence)
            else:
                yield name.split()[0], ''.join(sequence)


def get_fasta_size(filename):
    total_size = 0
    for _, seq in iterate_fasta(filename):
        total_size += len(seq)
    return total_size


def get_default_thread_count():
    return min(multiprocessing.cpu_count(), 16)


def get_n50(lengths):
    lengths = sorted(lengths, reverse=True)
    total_bases = sum(lengths)
    target_bases = total_bases * 0.5
    total_so_far = 0
    for length in lengths:
        total_so_far += length
        if total_so_far >= target_bases:
            return length
    return 0


def get_window_count(full_length, window_size, window_step):
    """
    Returns the number of possible windows in a sequence, given the full sequence length, the
    window size and the window step.

    For example:
      inputs: full_length=50, window_size=10, window_step=8
      output: 5 windows
    ----------      ----------      ----------
            ----------      ----------
    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    """
    if full_length < window_size:
        return 0
    return 1 + ((full_length - window_size) // window_step)


def get_window_coverage(window_size, window_step, window_count):
    """
    Returns the number of bases covered by sliding windows.

    For example:
      inputs: window_size=10, window_step=8, window_count=5
      output: 42
    ----------      ----------      ----------
            ----------      ----------
    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    """
    return window_size + ((window_count-1) * window_step)


def check_file_exists(filename):
    if filename.is_dir():
        sys.exit(f'Error: {filename} is a directory, not a file')
    if not filename.is_file():
        sys.exit(f'Error: {filename} does not exist')


def split_list(a, n):
    """
    https://stackoverflow.com/questions/2130016
    """
    k, m = divmod(len(a), n)
    return [a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n)]


def contains_ambiguous_bases(seq):
    unambiguous_bases = {'A', 'C', 'G', 'T'}
    return not all(base in unambiguous_bases for base in seq.upper())


def list_differences(a, b):
    """
    Assumes no duplicates in the lists.
    """
    a = set(a)
    b = set(b)
    in_both = a & b
    in_a_not_b = a - b
    in_b_not_a = b - a
    assert len(a | b) == len(in_both) + len(in_a_not_b) + len(in_b_not_a)
    return sorted(in_both), sorted(in_a_not_b), sorted(in_b_not_a)


def get_difference_count(cigar):
    """
    Returns the number of mismatches and indels in the CIGAR.
    """
    return cigar.count('X') + cigar.count('I') + cigar.count('D')
