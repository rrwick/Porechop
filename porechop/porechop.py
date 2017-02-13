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
import multiprocessing
from multiprocessing.dummy import Pool as ThreadPool
from .misc import load_fasta_or_fastq, print_table, red, bold_underline, check_file_exists, \
    MyHelpFormatter
from .adapters import ADAPTERS
from .nanopore_read import NanoporeRead
from .version import __version__


def main():
    args = get_arguments()
    reads, read_type = load_reads(args.input, args.verbosity, args.print_dest)

    matching_sets = find_matching_adapter_sets(reads, args.verbosity, args.end_size,
                                               args.scoring_scheme_vals, args.print_dest,
                                               args.adapter_threshold, args.check_reads,
                                               args.threads)

    display_adapter_set_results(matching_sets, args.verbosity, args.print_dest)

    if not matching_sets:
        print('No strong adapter matches - not proceeding with trimming\n', file=args.print_dest)
        sys.exit()

    find_adapters_at_read_ends(reads, matching_sets, args.verbosity, args.end_size,
                               args.extra_end_trim, args.end_threshold, args.scoring_scheme_vals,
                               args.print_dest, args.min_trim_size, args.threads)

    display_read_end_trimming_summary(reads, args.verbosity, args.print_dest)

    find_adapters_in_read_middles(reads, matching_sets, args.verbosity, args.middle_threshold,
                                  args.extra_middle_trim_good_side, args.extra_middle_trim_bad_side,
                                  args.scoring_scheme_vals, args.print_dest, args.threads)

    output_reads(reads, args.format, args.output, read_type, args.verbosity,
                 args.min_split_read_size, args.print_dest)


def get_arguments():
    """
    Parse the command line arguments.
    """
    default_threads = min(multiprocessing.cpu_count(), 16)

    parser = argparse.ArgumentParser(description='Porechop: a tool for finding adapters in Oxford '
                                                 'Nanopore reads, trimming them from the ends and '
                                                 'splitting reads with internal adapters',
                                     formatter_class=MyHelpFormatter)
    main_group = parser.add_argument_group('Main options')
    main_group.add_argument('-i', '--input', required=True,
                            help='FASTA or FASTQ of input reads (required)')
    main_group.add_argument('-o', '--output',
                            help='Filename for FASTA or FASTQ of trimmed reads (if not set, '
                                 'trimmed reads will be printed to stdout)')
    main_group.add_argument('--format', choices=['auto', 'fasta', 'fastq'], default='auto',
                            help='Output format for the reads - if auto, the '
                                 'format will be chosen based on the output filename or the input '
                                 'read format')
    main_group.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Level of progress information: 0 = none, 1 = some, 2 = full '
                                 ' - output will go to stdout if reads are saved to '
                                 'a file and stderr if reads are printed to stdout')
    main_group.add_argument('-t', '--threads', type=int, default=default_threads,
                            help='Number of threads to use for adapter alignment')

    main_group.add_argument('--version', action='version', version=__version__)

    adapter_search_group = parser.add_argument_group('Adapter search settings',
                                                     'Control how the program determines which '
                                                     'adapter sets are present')
    adapter_search_group.add_argument('--adapter_threshold', type=float, default=0.8,
                                      help='An adapter set has to score at least this well to be '
                                           'labelled as present and trimmed off (0.0 to 1.0)')
    adapter_search_group.add_argument('--check_reads', type=int, default=1000,
                                      help='This many reads will be aligned to all possible '
                                           'adapters to determine which adapter sets are present')
    adapter_search_group.add_argument('--scoring_scheme', type=str, default='3,-6,-5,-2',
                                      help='Comma-delimited string of alignment scores: match,'
                                           'mismatch, gap open, gap extend')

    end_trim_group = parser.add_argument_group('End adapter settings',
                                               'Control the trimming of adapters from read ends')
    end_trim_group.add_argument('--end_size', type=int, default=100,
                                help='The number of base pairs at each end of the read which will '
                                     'be searched for adapter sequences')
    end_trim_group.add_argument('--min_trim_size', type=int, default=4,
                                help='Adapter alignments smaller than this will be ignored')
    end_trim_group.add_argument('--extra_end_trim', type=int, default=2,
                                help='This many additional bases will be removed next to adapters '
                                     'found at the ends of reads')
    end_trim_group.add_argument('--end_threshold', type=float, default=0.5,
                                help='Adapters at the ends of reads must score at least this well '
                                     'to be removed (0.0 to 1.0)')

    middle_trim_group = parser.add_argument_group('Middle adapter settings',
                                                  'Control the splitting of read from middle '
                                                  'adapters')
    middle_trim_group.add_argument('--middle_threshold', type=float, default=0.7,
                                   help='Adapters in the middle of reads must score at least this '
                                        'well to be removed and split the read (0.0 to 1.0)')
    middle_trim_group.add_argument('--extra_middle_trim_good_side', type=int, default=10,
                                   help='This many additional bases will be removed next to '
                                        'middle adapters on their "good" side')
    middle_trim_group.add_argument('--extra_middle_trim_bad_side', type=int, default=100,
                                   help='This many additional bases will be removed next to '
                                        'middle adapters on their "bad" side')
    middle_trim_group.add_argument('--min_split_read_size', type=int, default=1000,
                                   help='Post-split read pieces smaller than this many base pairs '
                                        'will not be outputted')

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

    if args.threads < 1:
        sys.exit('Error: at least one thread required')

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


def find_matching_adapter_sets(reads, verbosity, end_size, scoring_scheme_vals, print_dest,
                               adapter_threshold, check_reads, threads):
    """
    Aligns all of the adapter sets to the start/end of reads to see which (if any) matches best.
    """
    if verbosity > 0:
        print(bold_underline('Looking for known adapter sets'), flush=True, file=print_dest)

    read_subset = reads[:check_reads]

    # If single-threaded, do the work in a simple loop.
    if threads == 1:
        for read in read_subset:
            for adapter_set in ADAPTERS:
                read.align_adapter_set(adapter_set, end_size, scoring_scheme_vals)

    # If multi-threaded, use a thread pool.
    else:
        def align_adapter_set_one_arg(all_args):
            r, a, b, c = all_args
            r.align_adapter_set(a, b, c)
        with ThreadPool(threads) as pool:
            arg_list = []
            for read in read_subset:
                for adapter_set in ADAPTERS:
                    arg_list.append((read, adapter_set, end_size, scoring_scheme_vals))
            pool.map(align_adapter_set_one_arg, arg_list)
    return [x for x in ADAPTERS if x.best_start_or_end_score() >= adapter_threshold]


def display_adapter_set_results(matching_sets, verbosity, print_dest):
    if verbosity < 1:
        return
    table = [['Set', 'Score']]
    row_colours = {}
    matching_set_names = [x.name for x in matching_sets]
    for adapter_set in ADAPTERS:
        score = adapter_set.best_start_or_end_score() * 100.0
        score = '%.1f' % score + '%'
        table.append([adapter_set.name, score])
        if adapter_set.name in matching_set_names:
            row_colours[len(table) - 1] = 'green'
    if verbosity > 0:
        print_table(table, print_dest, alignments='LR', row_colour=row_colours)
        print('\n', file=print_dest)


def find_adapters_at_read_ends(reads, matching_sets, verbosity, end_size, extra_trim_size,
                               end_threshold, scoring_scheme_vals, print_dest, min_trim_size,
                               threads):
    if verbosity > 0:
        matching_set_names = ', '.join([x.name for x in matching_sets])
        print(bold_underline('Trimming ' + matching_set_names + ' adapters from read ends'),
              file=print_dest)
        name_len = max(max(len(x.start_sequence[0]) for x in matching_sets),
                       max(len(x.end_sequence[0]) for x in matching_sets))
        for matching_set in matching_sets:
            print('  ' + matching_set.start_sequence[0].rjust(name_len) + ': ' +
                  red(matching_set.start_sequence[1]), file=print_dest)
            print('  ' + matching_set.end_sequence[0].rjust(name_len) + ': ' +
                  red(matching_set.end_sequence[1]), file=print_dest)
        print('', file=print_dest)

    # If single-threaded, do the work in a simple loop.
    if threads == 1:
        for read in reads:
            read.find_start_trim(matching_sets, end_size, extra_trim_size, end_threshold,
                                 scoring_scheme_vals, min_trim_size)
            read.find_end_trim(matching_sets, end_size, extra_trim_size, end_threshold,
                               scoring_scheme_vals, min_trim_size)
            if verbosity > 1:
                print(read.formatted_start_and_end_seq(end_size, extra_trim_size), file=print_dest)

    # If multi-threaded, use a thread pool.
    else:
        def start_end_trim_one_arg(all_args):
            r, a, b, c, d, e, f, v = all_args
            r.find_start_trim(a, b, c, d, e, f)
            r.find_end_trim(a, b, c, d, e, f)
            if v > 1:
                return r.formatted_start_and_end_seq(b, c)
            else:
                return ''
        with ThreadPool(threads) as pool:
            arg_list = []
            for read in reads:
                arg_list.append((read, matching_sets, end_size, extra_trim_size, end_threshold,
                                 scoring_scheme_vals, min_trim_size, verbosity))
            for out in pool.imap(start_end_trim_one_arg, arg_list):
                if verbosity > 1:
                    print(out, file=print_dest, flush=True)

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


def find_adapters_in_read_middles(reads, matching_sets, verbosity, middle_threshold,
                                  extra_trim_good_side, extra_trim_bad_side, scoring_scheme_vals,
                                  print_dest, threads):
    if verbosity > 0:
        matching_set_names = ', '.join([x.name for x in matching_sets])
        print(bold_underline('Splitting reads containing ' + matching_set_names + ' adapters'),
              file=print_dest)

    adapters = []
    for matching_set in matching_sets:
        adapters.append(matching_set.start_sequence)
        if matching_set.end_sequence[1] != matching_set.start_sequence[1]:
            adapters.append(matching_set.end_sequence)

    start_sequence_names = set(x.start_sequence[0] for x in matching_sets)
    end_sequence_names = set(x.end_sequence[0] for x in matching_sets)

    # If single-threaded, do the work in a simple loop.
    if threads == 1:
        for read in reads:
            read.find_middle_adapters(adapters, middle_threshold, extra_trim_good_side,
                                      extra_trim_bad_side, scoring_scheme_vals,
                                      start_sequence_names, end_sequence_names)

            if read.middle_adapter_positions and verbosity > 0:
                print(read.middle_adapter_results(verbosity), file=print_dest, flush=True)

    # If multi-threaded, use a thread pool.
    else:
        def find_middle_adapters_one_arg(all_args):
            r, a, b, c, d, e, f, g, v = all_args
            r.find_middle_adapters(a, b, c, d, e, f, g)
            return r.middle_adapter_results(v)
        with ThreadPool(threads) as pool:
            arg_list = []
            for read in reads:
                arg_list.append((read, adapters, middle_threshold, extra_trim_good_side,
                                 extra_trim_bad_side, scoring_scheme_vals, start_sequence_names,
                                 end_sequence_names, verbosity))
            for out in pool.imap(find_middle_adapters_one_arg, arg_list):
                if verbosity > 0 and out:
                    print(out, file=print_dest, flush=True)


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
            print(read_str, end='')

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
