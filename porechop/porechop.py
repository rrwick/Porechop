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

from __future__ import print_function, absolute_import
import argparse
import os
import sys
import subprocess
import multiprocessing
from multiprocessing.dummy import Pool as ThreadPool
from collections import defaultdict
from porechop.misc import load_fasta_or_fastq, print_table, red, bold_underline, check_file_exists, \
    MyHelpFormatter, int_to_str
from porechop.adapters import ADAPTERS
from porechop.nanopore_read import NanoporeRead
from porechop.version import __version__


def main():
    args = get_arguments()
    reads, read_type = load_reads(args.input, args.verbosity, args.print_dest)

    matching_sets = find_matching_adapter_sets(reads, args.verbosity, args.end_size,
                                               args.scoring_scheme_vals, args.print_dest,
                                               args.adapter_threshold, args.check_reads,
                                               args.threads)

    display_adapter_set_results(matching_sets, args.verbosity, args.print_dest)

    if matching_sets:
        check_barcodes = (args.barcode_dir is not None)
        find_adapters_at_read_ends(reads, matching_sets, args.verbosity, args.end_size,
                                   args.extra_end_trim, args.end_threshold,
                                   args.scoring_scheme_vals, args.print_dest, args.min_trim_size,
                                   args.threads, check_barcodes, args.barcode_threshold,
                                   args.barcode_diff, args.require_two_barcodes)
        display_read_end_trimming_summary(reads, args.verbosity, args.print_dest)

        find_adapters_in_read_middles(reads, matching_sets, args.verbosity, args.middle_threshold,
                                      args.extra_middle_trim_good_side,
                                      args.extra_middle_trim_bad_side, args.scoring_scheme_vals,
                                      args.print_dest, args.threads, args.discard_middle)
        display_read_middle_trimming_summary(reads, args.discard_middle, args.verbosity,
                                             args.print_dest)
    elif args.verbosity > 0:
        print('No adapters found - output reads are unchanged from input reads\n',
              file=args.print_dest)

    output_reads(reads, args.format, args.output, read_type, args.verbosity,
                 args.discard_middle, args.min_split_read_size, args.print_dest,
                 args.barcode_dir)


def get_arguments():
    """
    Parse the command line arguments.
    """
    default_threads = min(multiprocessing.cpu_count(), 16)

    parser = argparse.ArgumentParser(description=('Porechop: a tool for finding adapters in Oxford '
                                                 'Nanopore reads, trimming them from the ends and '
                                                 'splitting reads with internal adapters'),
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

    barcode_group = parser.add_argument_group('Barcode binning settings',
                                              'Control the binning of reads based on barcodes '
                                              '(i.e. barcode demultiplexing)')
    barcode_group.add_argument('-b', '--barcode_dir',
                               help='Reads will be binned based on their barcode and saved to '
                                    'separate files in this directory (incompatible with '
                                    '--output)')
    barcode_group.add_argument('--barcode_threshold', type=float, default=75.0,
                               help='A read must have at least this percent identity to a barcode '
                                    'to be binned')
    barcode_group.add_argument('--barcode_diff', type=float, default=5.0,
                               help="If the difference between a read's best barcode identity and "
                                    "its second-best barcode identity is less than this value, it "
                                    "will not be put in a barcode bin (to exclude cases which are "
                                    "too close to call)")
    barcode_group.add_argument('--require_two_barcodes', action='store_true',
                               help='Reads will only be put in barcode bins if they have a strong '
                                    'hit for the barcode on both their start and end (default: '
                                    'a read can be binned with only a single barcode alignment, '
                                    'assuming no contradictory barcode alignments exist at the '
                                    'other end)')

    adapter_search_group = parser.add_argument_group('Adapter search settings',
                                                     'Control how the program determines which '
                                                     'adapter sets are present')
    adapter_search_group.add_argument('--adapter_threshold', type=float, default=90.0,
                                      help='An adapter set has to have at least this percent '
                                           'identity to be labelled as present and trimmed off '
                                           '(0 to 100)')
    adapter_search_group.add_argument('--check_reads', type=int, default=10000,
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
    end_trim_group.add_argument('--end_threshold', type=float, default=75.0,
                                help='Adapters at the ends of reads must have at least this '
                                     'percent identity to be removed (0 to 100)')

    middle_trim_group = parser.add_argument_group('Middle adapter settings',
                                                  'Control the splitting of read from middle '
                                                  'adapters')
    middle_trim_group.add_argument('--discard_middle', action='store_true',
                                   help='Reads with middle adapters will be discarded (default: '
                                        'reads with middle adapters are split) (this option is '
                                        'on by default when outputting reads into barcode bins)')
    middle_trim_group.add_argument('--middle_threshold', type=float, default=85.0,
                                   help='Adapters in the middle of reads must have at least this '
                                        'percent identity to be found (0 to 100)')
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

    if args.barcode_dir is not None and args.output is not None:
        sys.exit('Error: only one of the following options may be used: --output, --barcode_dir')

    if args.barcode_dir is not None:
        args.discard_middle = True

    if args.output is None and args.barcode_dir is None:
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
    table = [['Set', 'Best %ID']]
    row_colours = {}
    matching_set_names = [x.name for x in matching_sets]
    for adapter_set in ADAPTERS:
        score = adapter_set.best_start_or_end_score()
        score = '%.1f' % score
        table.append([adapter_set.name, score])
        if adapter_set.name in matching_set_names:
            row_colours[len(table) - 1] = 'green'
    if verbosity > 0:
        print_table(table, print_dest, alignments='LR', row_colour=row_colours)
        print('\n', file=print_dest)


def find_adapters_at_read_ends(reads, matching_sets, verbosity, end_size, extra_trim_size,
                               end_threshold, scoring_scheme_vals, print_dest, min_trim_size,
                               threads, check_barcodes, barcode_threshold, barcode_diff,
                               require_two_barcodes):
    if verbosity > 0:
        print(bold_underline('Trimming adapters from read ends'),
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
                                 scoring_scheme_vals, min_trim_size, check_barcodes)
            read.find_end_trim(matching_sets, end_size, extra_trim_size, end_threshold,
                               scoring_scheme_vals, min_trim_size, check_barcodes)
            if check_barcodes:
                read.determine_barcode(barcode_threshold, barcode_diff, require_two_barcodes)
            if verbosity > 1:
                print(read.formatted_start_and_end_seq(end_size, extra_trim_size, check_barcodes),
                      file=print_dest)

    # If multi-threaded, use a thread pool.
    else:
        def start_end_trim_one_arg(all_args):
            r, a, b, c, d, e, f, g, h, i, j, v = all_args
            r.find_start_trim(a, b, c, d, e, f, g)
            r.find_end_trim(a, b, c, d, e, f, g)
            if check_barcodes:
                r.determine_barcode(h, i, j)
            if v > 1:
                return r.formatted_start_and_end_seq(b, c, g)
            else:
                return ''
        with ThreadPool(threads) as pool:
            arg_list = []
            for read in reads:
                arg_list.append((read, matching_sets, end_size, extra_trim_size, end_threshold,
                                 scoring_scheme_vals, min_trim_size, check_barcodes,
                                 barcode_threshold, barcode_diff, require_two_barcodes, verbosity))
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
                                  print_dest, threads, discard_middle):
    if verbosity > 0:
        verb = 'Discarding' if discard_middle else 'Splitting'
        print(bold_underline(verb + ' reads containing middle adapters'),
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


def display_read_middle_trimming_summary(reads, discard_middle, verbosity, print_dest):
    if verbosity < 1:
        return
    middle_trim_count = sum(1 if x.middle_adapter_positions else 0 for x in reads)
    verb = 'discarded' if discard_middle else 'split'
    print(str(middle_trim_count) + ' / ' + str(len(reads)) + ' reads were ' + verb +
          ' based on middle adapters\n', file=print_dest)


def output_reads(reads, out_format, output, read_type, verbosity, discard_middle,
                 min_split_size, print_dest, barcode_dir):
    if out_format == 'auto':
        if output is None:
            out_format = read_type.lower()
        elif '.fasta' in output:
            out_format = 'fasta'
        elif '.fastq' in output:
            out_format = 'fastq'
        else:
            out_format = read_type.lower()

    # Output reads to barcode bins.
    if barcode_dir is not None:
        if not os.path.isdir(barcode_dir):
            os.makedirs(barcode_dir)
        barcode_files = {}
        barcode_read_counts, barcode_base_counts = defaultdict(int), defaultdict(int)
        for read in reads:
            barcode_name = read.barcode_call
            read_str = read.get_fasta(min_split_size, discard_middle) if out_format == 'fasta' \
                else read.get_fastq(min_split_size, discard_middle)
            if not read_str:
                continue
            if barcode_name not in barcode_files:
                barcode_files[barcode_name] = \
                    open(os.path.join(barcode_dir, barcode_name + '.' + out_format), 'wt')
            barcode_files[barcode_name].write(read_str)
            barcode_read_counts[barcode_name] += 1
            barcode_base_counts[barcode_name] += read.seq_length_with_start_end_adapters_trimmed()
        table = [['Barcode', 'Reads', 'Bases', 'File']]

        for barcode_name in sorted(barcode_files.keys()):
            barcode_files[barcode_name].close()
            bin_filename = os.path.join(barcode_dir, barcode_name + '.' + out_format)
            if not os.path.isfile(bin_filename):
                continue
            bin_filename_gz = bin_filename + '.gz'
            if os.path.isfile(bin_filename_gz):
                os.remove(bin_filename_gz)
            subprocess.check_output('gzip ' + bin_filename, stderr=subprocess.STDOUT, shell=True)
            table_row = [barcode_name, int_to_str(barcode_read_counts[barcode_name]),
                         int_to_str(barcode_base_counts[barcode_name]), bin_filename_gz]
            table.append(table_row)
        if verbosity > 0:
            print('')
            print(bold_underline('Saving trimmed reads to barcode-specific files'), flush=True,
                  file=print_dest)
            print_table(table, print_dest, alignments='LRRL', max_col_width=60, col_separation=2)
            print('')

    # Output to all reads to stdout.
    elif output is None:
        for read in reads:
            read_str = read.get_fasta(min_split_size, discard_middle) if out_format == 'fasta' \
                else read.get_fastq(min_split_size, discard_middle)
            print(read_str, end='')

    # Output to all reads to file.
    else:
        gzipped_out = output.endswith('.gz')
        if gzipped_out:
            out_filename = 'TEMP_' + str(os.getpid()) + '.fastq'
        else:
            out_filename = output
        with open(out_filename, 'wt') as out:
            for read in reads:
                read_str = read.get_fasta(min_split_size, discard_middle) if out_format == 'fasta' \
                    else read.get_fastq(min_split_size, discard_middle)
                out.write(read_str)
        if gzipped_out:
            subprocess.check_output('gzip -c ' + out_filename + ' > ' + output,
                                    stderr=subprocess.STDOUT, shell=True)
            os.remove(out_filename)
        if verbosity > 0:
            print('\nSaved result to ' + os.path.abspath(output) + '\n', file=print_dest)
