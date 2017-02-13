#!/usr/bin/env python3
"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Porechop

This module contains the main script for Porechop. It is executed when a user runs `porechop`
(after installation) or `porechop-runner.py` (directly from the source directory).

This file is part of Porechop. Porechop is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Porechop is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Porechop. If
not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import os
import sys
import subprocess
from .misc import load_fasta_or_fastq, print_table, red, bold_underline, check_file_exists
from .adapters import ADAPTERS
from .nanopore_read import NanoporeRead
from .version import __version__


def main():
    args = get_arguments()
    reads, read_type = load_reads(args.input, args.verbosity, args.print_dest)

    best_adapter = find_best_match_adapter_set(reads, args.verbosity, args.end_size,
                                               args.scoring_scheme_vals, args.print_dest)

    display_adapter_set_results(best_adapter, args.verbosity, args.adapter_threshold,
                                args.print_dest)

    if best_adapter.get_average_best_score() < args.adapter_threshold:
        print('No strong adapter matches - not proceeding with trimming\n', file=args.print_dest)
        sys.exit()

    find_adapters_at_read_ends(reads, best_adapter, args.verbosity, args.end_size,
                               args.extra_end_trim, args.end_threshold, args.scoring_scheme_vals,
                               args.print_dest)

    display_read_end_trimming_summary(reads, args.verbosity, args.print_dest)

    find_adapters_in_read_middles(reads, best_adapter, args.verbosity, args.middle_threshold,
                                  args.extra_middle_trim_good_side, args.extra_middle_trim_bad_side,
                                  args.scoring_scheme_vals, args.print_dest)

    output_reads(reads, args.format, args.output, read_type, args.verbosity,
                 args.min_split_read_size, args.print_dest)


def get_arguments():
    """
    Parse the command line arguments.
    """
    parser = argparse.ArgumentParser(description='Porechop: a tool for finding adapters in Oxford '
                                                 'Nanopore reads, trimming them from the ends and '
                                                 'splitting reads with internal adapters')
    main_group = parser.add_argument_group('Main options')
    main_group.add_argument('-i', '--input', required=True,
                            help='FASTA or FASTQ of input reads (required)')
    main_group.add_argument('-o', '--output',
                            help='Filename for FASTA or FASTQ of trimmed reads (if not set, '
                                 'trimmed reads will be printed to stdout)')
    main_group.add_argument('--format', choices=['auto', 'fasta', 'fastq'], default='auto',
                            help='Output format for the reads (default: auto) - if auto, the '
                                 'format will be chosen based on the output filename or the input '
                                 'read format)')
    main_group.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Level of progress information: 0 = none, 1 = some, 2 = full '
                                 '(default: 1) - info will be sent to stdout if reads are saved to '
                                 'a file and stderr if reads are printed to stdout')
    main_group.add_argument('--version', action='version', version=__version__)

    adapter_search_group = parser.add_argument_group('Adapter search settings',
                                                     'Control how the program determines which '
                                                     'adapter sets are present')
    adapter_search_group.add_argument('--adapter_threshold', type=float, default=0.7,
                                      help='An adapter set has to score at least this well to be '
                                           'labelled as present and trimmed off (0.0 to 1.0, '
                                           'default: 0.7)')
    adapter_search_group.add_argument('--scoring_scheme', type=str, default='3,-6,-5,-2',
                                      help='Comma-delimited string of alignment scores: match,'
                                           'mismatch, gap open, gap extend (default: "3,-6,-5,-2"')

    end_trim_group = parser.add_argument_group('End adapter settings',
                                               'Control the trimming of adapters from read ends')
    end_trim_group.add_argument('--end_size', type=int, default=60,
                                help='The number of base pairs at each end of the read which will '
                                     'be searched for adapter sequences')
    end_trim_group.add_argument('--extra_end_trim', type=int, default=2,
                                help='This many additional bases will be removed next to adapters '
                                     'found at the ends of reads')
    end_trim_group.add_argument('--end_threshold', type=float, default=0.5,
                                help='Adapters at the ends of reads must score at least this well '
                                     'to be removed (0.0 to 1.0, default: 0.5)')

    middle_trim_group = parser.add_argument_group('Middle adapter settings',
                                                  'Control the splitting of read from middle '
                                                  'adapters')
    middle_trim_group.add_argument('--middle_threshold', type=float, default=0.7,
                                   help='Adapters in the middle of reads must score at least this '
                                        'well to be removed and split the read (0.0 to 1.0, '
                                        'default: 0.7)')
    middle_trim_group.add_argument('--extra_middle_trim_good_side', type=int, default=10,
                                   help='This many additional bases will be removed next to '
                                        'middle adapters on their "good" side (default: 10)')
    middle_trim_group.add_argument('--extra_middle_trim_bad_side', type=int, default=100,
                                   help='This many additional bases will be removed next to '
                                        'middle adapters on their "bad" side (default: 100)')
    middle_trim_group.add_argument('--min_split_read_size', type=int, default=1000,
                                   help='Post-split read pieces smaller than this many base pairs '
                                        'will not be outputted (default: 1000)')

    args = parser.parse_args()

    try:
        scoring_scheme = [int(x) for x in args.scoring_scheme.split(',')]
    except ValueError:
        sys.exit('Error: incorrectly formatted scoring scheme')
    if len(scoring_scheme) != 4:
        sys.exit('Error: incorrectly formatted scoring scheme')
    args.scoring_scheme_vals = scoring_scheme

    if args.output is None:
        args.print_dest = sys.stderr
    else:
        args.print_dest = sys.stdout

    return args


def load_reads(input_filename, verbosity, print_dest):
    check_file_exists(input_filename)
    reads, read_type = load_fasta_or_fastq(input_filename)
    if read_type == 'FASTA':
        reads = [NanoporeRead(x[2], x[1], '') for x in reads]
    else:  # FASTQ
        reads = [NanoporeRead(x[4], x[1], x[3]) for x in reads]
    if verbosity > 0:
        print('', file=print_dest)
    return reads, read_type


def find_best_match_adapter_set(reads, verbosity, end_size, scoring_scheme_vals, print_dest):
    """
    Aligns all of the adapter sets to the start/end of reads to see which (if any) matches best.
    """
    if verbosity > 0:
        print(bold_underline('Looking for known adapter sets'), flush=True, file=print_dest)
    for read in reads:
        for adapter_set in ADAPTERS:
            read.align_adapter_set(adapter_set, end_size, scoring_scheme_vals)
    best_adapter = sorted(ADAPTERS, key=lambda x: x.get_average_best_score())[-1]
    return best_adapter


def display_adapter_set_results(best_adapter, verbosity, adapter_threshold, print_dest):
    if verbosity < 1:
        return
    table = [['Adapter set', 'Score']]
    best_row_num = -1
    for adapter_set in ADAPTERS:
        score = adapter_set.get_average_best_score() * 100.0
        score = '%.1f' % score + '%'
        table.append([adapter_set.name, score])
        if adapter_set.name == best_adapter.name and \
                adapter_set.get_average_best_score() >= adapter_threshold:
            best_row_num = len(table) - 1
    if verbosity > 0:
        print_table(table, print_dest, alignments='LR', row_colour={best_row_num: 'green'})
        print('\n', file=print_dest)


def find_adapters_at_read_ends(reads, best_adapter, verbosity, end_size, extra_trim_size,
                               end_threshold, scoring_scheme_vals, print_dest):
    if verbosity > 0:
        print(bold_underline('Trimming ' + best_adapter.name + ' adapters from read ends'),
              file=print_dest)
        name_len = max(len(best_adapter.start_sequence[0]), len(best_adapter.end_sequence[0]))
        print('  ' + best_adapter.start_sequence[0].rjust(name_len) + ': ' +
              red(best_adapter.start_sequence[1]), file=print_dest)
        print('  ' + best_adapter.end_sequence[0].rjust(name_len) + ': ' +
              red(best_adapter.end_sequence[1]) + '\n', file=print_dest)
    for read in reads:
        read.find_start_trim(best_adapter, end_size, extra_trim_size, end_threshold,
                             scoring_scheme_vals)
        read.find_end_trim(best_adapter, end_size, extra_trim_size, end_threshold,
                           scoring_scheme_vals)
        if verbosity > 1:
            print(read.get_formatted_start_seq(end_size, extra_trim_size) + '...' +
                  read.get_formatted_end_seq(end_size, extra_trim_size), file=print_dest)
    if verbosity > 1:
        print('', file=print_dest)


def display_read_end_trimming_summary(reads, verbosity, print_dest):
    if verbosity < 1:
        return
    start_trim_total = sum(x.start_trim_amount for x in reads)
    start_trim_count = sum(1 if x.start_trim_amount else 0 for x in reads)
    end_trim_count = sum(1 if x.end_trim_amount else 0 for x in reads)
    end_trim_total = sum(x.end_trim_amount for x in reads)
    print(str(start_trim_count).rjust(len(str(len(reads)))) + ' / ' + str(len(reads)) +
          ' reads had adapters trimmed from their start (' + str(start_trim_total) +
          ' bp removed)', file=print_dest)
    print(str(end_trim_count).rjust(len(str(len(reads)))) + ' / ' + str(len(reads)) +
          ' reads had adapters trimmed from their end (' + str(end_trim_total) +
          ' bp removed)', file=print_dest)
    print('\n', file=print_dest)


def find_adapters_in_read_middles(reads, best_adapter, verbosity, middle_threshold,
                                  extra_middle_trim_good_side, extra_middle_trim_bad_side,
                                  scoring_scheme_vals, print_dest):
    print(bold_underline('Splitting reads containing ' + best_adapter.name + ' adapters'),
          file=print_dest)
    for read in reads:
        read.find_middle_adapters(best_adapter, verbosity, middle_threshold,
                                  extra_middle_trim_good_side, extra_middle_trim_bad_side,
                                  scoring_scheme_vals, print_dest)


def output_reads(reads, out_format, output, read_type, verbosity, min_split_read_size, print_dest):
    if out_format == 'auto':
        if output is None:
            out_format = read_type
        elif '.fasta' in output:
            out_format = 'fasta'
        elif '.fastq' in output:
            out_format = 'fastq'
        else:
            out_format = read_type.lower()

    if output is None:  # output to stdout
        for read in reads:
            read_str = read.get_fasta(min_split_read_size) if out_format == 'fasta' \
                else read.get_fastq(min_split_read_size)
            print(read_str, end='', file=print_dest)

    else:  # output to file
        gzipped_out = output.endswith('.gz')
        if gzipped_out:
            out_filename = 'TEMP_' + str(os.getpid()) + '.fastq'
        else:
            out_filename = output
        with open(out_filename, 'wt') as out:
            for read in reads:
                read_str = read.get_fasta() if out_format == 'fasta' else read.get_fastq()
                out.write(read_str)
        if gzipped_out:
            subprocess.check_output('gzip -c ' + out_filename + ' > ' + output,
                                    stderr=subprocess.STDOUT, shell=True)
            os.remove(out_filename)
        if verbosity > 1:
            print('\nSaved result to ' + os.path.abspath(output) + '\n', file=print_dest)
