#!/usr/bin/env python3
"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Porechop

This module contains the class for a Nanopore read.

This file is part of Porechop. Porechop is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Porechop is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Porechop. If
not, see <http://www.gnu.org/licenses/>.
"""

from .cpp_function_wrappers import adapter_alignment
from .misc import yellow, red, add_line_breaks_to_sequence, END_FORMATTING, RED, YELLOW


class NanoporeRead(object):

    def __init__(self, name, seq, quals):
        self.name = name

        self.seq = seq
        self.quals = quals
        if not self.quals:
            self.quals = '+' * len(seq)

        self.start_trim_amount = 0
        self.end_trim_amount = 0

        self.middle_adapter_positions = set()
        self.middle_trim_positions = set()

    def get_seq_with_start_end_adapters_trimmed(self):
        if not self.start_trim_amount and not self.end_trim_amount:
            return self.seq
        start_pos = self.start_trim_amount
        end_pos = len(self.seq) - self.end_trim_amount
        trimmed_seq = self.seq[start_pos:end_pos]
        return trimmed_seq

    def get_quals_with_start_end_adapters_trimmed(self):
        if not self.start_trim_amount and not self.end_trim_amount:
            return self.quals
        start_pos = self.start_trim_amount
        end_pos = len(self.quals) - self.end_trim_amount
        trimmed_quals = self.quals[start_pos:end_pos]
        return trimmed_quals

    def get_split_read_parts(self, min_split_read_size):
        """
        Returns the read split into parts as determined by the middle_trim_positions set.
        """
        trimmed_seq = self.get_seq_with_start_end_adapters_trimmed()
        trimmed_quals = self.get_quals_with_start_end_adapters_trimmed()
        split_read_parts = []
        part_seq, part_quals = [], []
        for i in range(len(trimmed_seq)):
            if i in self.middle_trim_positions:
                if part_seq:
                    split_read_parts.append((''.join(part_seq), ''.join(part_quals)))
                    part_seq, part_quals = [], []
            else:
                part_seq.append(trimmed_seq[i])
                part_quals.append(trimmed_quals[i])
        if part_seq:
            split_read_parts.append((''.join(part_seq), ''.join(part_quals)))
        split_read_parts = [x for x in split_read_parts if len(x[0]) >= min_split_read_size]
        return split_read_parts

    def get_fasta(self, min_split_read_size):
        if not self.middle_trim_positions:
            seq = add_line_breaks_to_sequence(self.get_seq_with_start_end_adapters_trimmed(), 70)
            return ''.join(['>', self.name, '\n', seq])
        else:
            fasta_str = ''
            for i, split_read_part in enumerate(self.get_split_read_parts(min_split_read_size)):
                read_name = add_number_to_read_name(self.name, i + 1)
                seq = add_line_breaks_to_sequence(split_read_part[0], 70)
                fasta_str += ''.join(['>', read_name, '\n', seq])
            return fasta_str

    def get_fastq(self, min_split_read_size):
        if not self.middle_trim_positions:
            return ''.join(['@', self.name, '\n',
                            self.get_seq_with_start_end_adapters_trimmed(), '\n+\n',
                            self.get_quals_with_start_end_adapters_trimmed(), '\n'])
        else:
            fastq_str = ''
            for i, split_read_part in enumerate(self.get_split_read_parts(min_split_read_size)):
                read_name = add_number_to_read_name(self.name, i + 1)
                fastq_str += ''.join(['@', read_name, '\n', split_read_part[0], '\n+\n',
                                      split_read_part[1], '\n'])
            return fastq_str

    def align_adapter_set(self, adapter_set, end_size, scoring_scheme_vals):
        """
        This function aligns the adapter to the reads and updates the best score for the adapter.
        This is not to determine where to trim the reads, but rather to figure out which adapter
        sets are present in the data.
        """
        read_seq_start = self.seq[:end_size]
        score, _, _, _ = align_adapter(read_seq_start, adapter_set.start_sequence[1],
                                       scoring_scheme_vals)
        adapter_set.best_start_score = max(adapter_set.best_start_score, score)
        read_seq_end = self.seq[-end_size:]
        score, _, _, _ = align_adapter(read_seq_end, adapter_set.end_sequence[1],
                                       scoring_scheme_vals)
        adapter_set.best_end_score = max(adapter_set.best_end_score, score)

    def find_start_trim(self, adapter_set, end_size, extra_trim_size, end_threshold,
                        scoring_scheme_vals, min_trim_size):
        """
        Aligns an adapter sequence and possibly adjusts the read's start trim amount based on the
        result.
        """
        read_seq_start = self.seq[:end_size]
        _, score, read_start, read_end = align_adapter(read_seq_start,
                                                       adapter_set.start_sequence[1],
                                                       scoring_scheme_vals)
        if score > end_threshold and read_end != end_size and \
                read_end - read_start >= min_trim_size:
            trim_amount = read_end + extra_trim_size
            self.start_trim_amount = max(self.start_trim_amount, trim_amount)

    def find_end_trim(self, adapter_set, end_size, extra_trim_size, end_threshold,
                      scoring_scheme_vals, min_trim_size):
        """
        Aligns an adapter sequence and possibly adjusts the read's end trim amount based on the
        result.
        """
        read_seq_end = self.seq[-end_size:]
        _, score, read_start, read_end = align_adapter(read_seq_end, adapter_set.end_sequence[1],
                                                       scoring_scheme_vals)
        if score > end_threshold and read_start != 0 and \
                read_end - read_start >= min_trim_size:
            trim_amount = (end_size - read_start) + extra_trim_size
            self.end_trim_amount = max(self.end_trim_amount, trim_amount)

    def find_middle_adapters(self, adapter_set, verbosity, middle_threshold,
                             extra_middle_trim_good_side, extra_middle_trim_bad_side,
                             scoring_scheme_vals, print_dest):
        """
        Aligns an adapter sequence to the whole read to find places where the read should be split.
        """
        adapters = [adapter_set.start_sequence]
        if adapter_set.end_sequence[1] != adapter_set.start_sequence[1]:
            adapters += [adapter_set.end_sequence]

        masked_seq = self.get_seq_with_start_end_adapters_trimmed()
        first_hit = True

        for adapter_name, adapter_seq in adapters:

            # We keep aligning adapters as long we get strong hits, so we can find multiple
            # occurrences in a single read.
            while True:
                score, _, read_start, read_end = align_adapter(masked_seq, adapter_seq,
                                                               scoring_scheme_vals)
                if score >= middle_threshold:
                    if first_hit:
                        print(self.name, file=print_dest)
                        first_hit = False

                    print('  found ' + adapter_name + ' (pos ' + str(read_start) + '-' +
                          str(read_end) + ')', file=print_dest)

                    masked_seq = masked_seq[:read_start] + '-' * (read_end - read_start) + \
                        masked_seq[read_end:]
                    self.middle_adapter_positions.update(range(read_start, read_end))

                    trim_start = read_start - extra_middle_trim_good_side
                    if adapter_name == adapter_set.start_sequence[0]:
                        trim_start = read_start - extra_middle_trim_bad_side
                    trim_end = read_end + extra_middle_trim_good_side
                    if adapter_name == adapter_set.end_sequence[0]:
                        trim_end = read_end + extra_middle_trim_bad_side

                    self.middle_trim_positions.update(range(trim_start, trim_end))
                else:
                    break

        if self.middle_adapter_positions:
            if verbosity > 1:
                print(self.get_formatted_middle_seq(), file=print_dest)
            print('', file=print_dest)

    def get_formatted_start_seq(self, end_size, extra_trim_size):
        """
        Returns the start of the read sequence, with any found adapters highlighted in red.
        """
        start_seq = self.seq[:end_size]
        if not self.start_trim_amount:
            return start_seq
        red_bases = self.start_trim_amount - extra_trim_size
        formatted_str = ''
        if red_bases:
            formatted_str = red(start_seq[:red_bases])
        formatted_str += yellow(start_seq[red_bases:red_bases+2])
        formatted_str += start_seq[red_bases+2:]
        return formatted_str

    def get_formatted_end_seq(self, end_size, extra_trim_size):
        """
        Returns the end of the read sequence, with any found adapters highlighted in red.
        """
        end_seq = self.seq[-end_size:]
        if not self.end_trim_amount:
            return end_seq
        red_bases = self.end_trim_amount - extra_trim_size
        formatted_str = ''
        if red_bases:
            formatted_str = red(end_seq[-red_bases:])
        formatted_str = yellow(end_seq[-(red_bases+2):-red_bases]) + formatted_str
        formatted_str = end_seq[:-(red_bases+2)] + formatted_str
        return formatted_str

    def get_formatted_middle_seq(self):
        """
        If a middle adapter was found, this returns the relevant part of the read sequence, with
        the adapter highlighted in red.
        """
        if not self.middle_adapter_positions:
            return

        trimmed_seq = self.get_seq_with_start_end_adapters_trimmed()

        range_start = max(0, min(self.middle_trim_positions) - 100)
        range_end = min(len(trimmed_seq),
                        max(self.middle_trim_positions) + 100)
        formatted_str = '' if range_start == 0 else '(' + str(range_start) + ' bp)...'

        last_colour = None
        for i in range(range_start, range_end):
            char_colour = None
            if i in self.middle_trim_positions:
                char_colour = 'yellow'
            if i in self.middle_adapter_positions:
                char_colour = 'red'
            if char_colour != last_colour:
                formatted_str += END_FORMATTING
                if char_colour == 'yellow':
                    formatted_str += YELLOW
                if char_colour == 'red':
                    formatted_str += RED

            formatted_str += trimmed_seq[i]
            last_colour = char_colour
        if last_colour is not None:
            formatted_str += END_FORMATTING

        formatted_str += '' if range_end == len(trimmed_seq) \
            else '...(' + str(len(trimmed_seq) - range_end) + ' bp)'
        return formatted_str


def align_adapter(read_seq, adapter_seq, scoring_scheme_vals):
    alignment_result = adapter_alignment(read_seq, adapter_seq, scoring_scheme_vals)

    result_parts = alignment_result.split(',')

    read_start = int(result_parts[0])
    read_end = int(result_parts[1]) + 1
    adapter_start = int(result_parts[2])
    adapter_end = int(result_parts[3]) + 1
    raw_score = int(result_parts[4])

    aligned_adapter_length = abs(adapter_end - adapter_start)

    try:
        full_length_score = raw_score / (scoring_scheme_vals[0] * len(adapter_seq))
    except ZeroDivisionError:
        full_length_score = 0
    try:
        aligned_length_score = raw_score / (scoring_scheme_vals[0] * aligned_adapter_length)
    except ZeroDivisionError:
        aligned_length_score = 0
    return full_length_score, aligned_length_score, read_start, read_end


def add_number_to_read_name(read_name, number):
    if ' ' not in read_name:
        return read_name + '_' + str(number)
    else:
        return read_name.replace(' ', '_' + str(number) + ' ', 1)
