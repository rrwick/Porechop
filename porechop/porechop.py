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
import shutil
import re
from multiprocessing.dummy import Pool as ThreadPool
from collections import defaultdict
from .misc import load_fasta_or_fastq, print_table, red, bold_underline, MyHelpFormatter, int_to_str
from .adapters import ADAPTERS, make_full_native_barcode_adapter,\
    make_old_full_rapid_barcode_adapter, make_new_full_rapid_barcode_adapter
from .nanopore_read import NanoporeRead
from .version import __version__


def main():
    args = get_arguments()
    reads, check_reads, read_type = load_reads(args.input, args.verbosity, args.print_dest,
                                               args.check_reads)

    matching_sets = find_matching_adapter_sets(check_reads, args.verbosity, args.end_size,
                                               args.scoring_scheme_vals, args.print_dest,
                                               args.adapter_threshold, args.threads)
    matching_sets = fix_up_1d2_sets(matching_sets)

    if args.barcode_dir:
        forward_or_reverse_barcodes = choose_barcoding_kit(matching_sets, args.verbosity,
                                                           args.print_dest)
    else:
        forward_or_reverse_barcodes = None

    display_adapter_set_results(matching_sets, args.verbosity, args.print_dest)
    matching_sets = add_full_barcode_adapter_sets(matching_sets)

    if args.verbosity > 0:
        print('\n', file=args.print_dest)

    if matching_sets:
        check_barcodes = (args.barcode_dir is not None)
        find_adapters_at_read_ends(reads, matching_sets, args.verbosity, args.end_size,
                                   args.extra_end_trim, args.end_threshold,
                                   args.scoring_scheme_vals, args.print_dest, args.min_trim_size,
                                   args.threads, check_barcodes, args.barcode_threshold,
                                   args.barcode_diff, args.require_two_barcodes,
                                   forward_or_reverse_barcodes)
        display_read_end_trimming_summary(reads, args.verbosity, args.print_dest)

        if not args.no_split:
            find_adapters_in_read_middles(reads, matching_sets, args.verbosity,
                                          args.middle_threshold, args.extra_middle_trim_good_side,
                                          args.extra_middle_trim_bad_side, args.scoring_scheme_vals,
                                          args.print_dest, args.threads, args.discard_middle)
            display_read_middle_trimming_summary(reads, args.discard_middle, args.verbosity,
                                                 args.print_dest)
    elif args.verbosity > 0:
        print('No adapters found - output reads are unchanged from input reads\n',
              file=args.print_dest)

    output_reads(reads, args.format, args.output, read_type, args.verbosity,
                 args.discard_middle, args.min_split_read_size, args.print_dest,
                 args.barcode_dir, args.input, args.untrimmed, args.threads,
                 args.discard_unassigned)


def get_arguments():
    """
    Parse the command line arguments.
    """
    default_threads = min(multiprocessing.cpu_count(), 16)

    parser = argparse.ArgumentParser(description='Porechop: a tool for finding adapters in Oxford '
                                                 'Nanopore reads, trimming them from the ends and '
                                                 'splitting reads with internal adapters',
                                     formatter_class=MyHelpFormatter, add_help=False)
    main_group = parser.add_argument_group('Main options')
    main_group.add_argument('-i', '--input', required=True,
                            help='FASTA/FASTQ of input reads or a directory which will be '
                                 'recursively searched for FASTQ files (required)')
    main_group.add_argument('-o', '--output',
                            help='Filename for FASTA or FASTQ of trimmed reads (if not set, '
                                 'trimmed reads will be printed to stdout)')
    main_group.add_argument('--format', choices=['auto', 'fasta', 'fastq', 'fasta.gz', 'fastq.gz'],
                            default='auto',
                            help='Output format for the reads - if auto, the '
                                 'format will be chosen based on the output filename or the input '
                                 'read format')
    main_group.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Level of progress information: 0 = none, 1 = some, 2 = lots, '
                                 '3 = full - output will go to stdout if reads are saved to '
                                 'a file and stderr if reads are printed to stdout')
    main_group.add_argument('-t', '--threads', type=int, default=default_threads,
                            help='Number of threads to use for adapter alignment')

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
                                    'match for the barcode on both their start and end (default: '
                                    'a read can be binned with a match at its start or end)')
    barcode_group.add_argument('--untrimmed', action='store_true',
                               help='Bin reads but do not trim them (default: trim the reads)')
    barcode_group.add_argument('--discard_unassigned', action='store_true',
                               help='Discard unassigned reads (instead of creating a "none" bin)')

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
                                      help='Comma-delimited string of alignment scores: match, '
                                           'mismatch, gap open, gap extend')

    end_trim_group = parser.add_argument_group('End adapter settings',
                                               'Control the trimming of adapters from read ends')
    end_trim_group.add_argument('--end_size', type=int, default=150,
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
    middle_trim_group.add_argument('--no_split', action='store_true',
                                   help='Skip splitting reads based on middle adapters '
                                        '(default: split reads when an adapter is found in the '
                                        'middle)')
    middle_trim_group.add_argument('--discard_middle', action='store_true',
                                   help='Reads with middle adapters will be discarded (default: '
                                        'reads with middle adapters are split) (required for '
                                        'reads to be used with Nanopolish, this option is on by '
                                        'default when outputting reads into barcode bins)')
    middle_trim_group.add_argument('--middle_threshold', type=float, default=90.0,
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

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version', version=__version__,
                           help="Show program's version number and exit")

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

    if args.untrimmed and args.barcode_dir is None:
        sys.exit('Error: --untrimmed can only be used with --barcode_dir')

    if args.barcode_dir is not None:
        args.discard_middle = True

    if args.output is None and args.barcode_dir is None:
        args.print_dest = sys.stderr
    else:
        args.print_dest = sys.stdout

    if args.threads < 1:
        sys.exit('Error: at least one thread required')

    return args


def load_reads(input_file_or_directory, verbosity, print_dest, check_read_count):

    # If the input is a file, just load reads from that file. The check reads will just be the
    # first reads from that file.
    if os.path.isfile(input_file_or_directory):
        if verbosity > 0:
            print('\n' + bold_underline('Loading reads'), flush=True, file=print_dest)
            print(input_file_or_directory, flush=True, file=print_dest)
        reads, read_type = load_fasta_or_fastq(input_file_or_directory)
        if read_type == 'FASTA':
            reads = [NanoporeRead(x[2], x[1], '') for x in reads]
        else:  # FASTQ
            reads = [NanoporeRead(x[4], x[1], x[3]) for x in reads]
        check_reads = reads[:check_read_count]

    # If the input is a directory, assume it's an Albacore directory and search it recursively for
    # fastq files. The check reads will be spread over all of the input files.
    elif os.path.isdir(input_file_or_directory):
        if verbosity > 0:
            print('\n' + bold_underline('Searching for FASTQ files'), flush=True, file=print_dest)
        fastqs = sorted([os.path.join(dir_path, f)
                         for dir_path, _, filenames in os.walk(input_file_or_directory)
                         for f in filenames
                         if f.lower().endswith('.fastq') or f.lower().endswith('.fastq.gz')])
        if not fastqs:
            sys.exit('Error: could not find fastq files in ' + input_file_or_directory)
        reads = []
        read_type = 'FASTQ'
        check_reads = []
        check_reads_per_file = int(round(check_read_count / len(fastqs)))
        for fastq_file in fastqs:
            if verbosity > 0:
                print(fastq_file, flush=True, file=print_dest)
            file_reads, _ = load_fasta_or_fastq(fastq_file)
            file_reads = [NanoporeRead(x[4], x[1], x[3]) for x in file_reads]

            albacore_barcode = get_albacore_barcode_from_path(fastq_file)
            for read in file_reads:
                read.albacore_barcode_call = albacore_barcode
            reads += file_reads
            check_reads += file_reads[:check_reads_per_file]
        if verbosity > 0:
            print('', flush=True, file=print_dest)

    else:
        sys.exit('Error: could not find ' + input_file_or_directory)

    if verbosity > 0:
        print(int_to_str(len(reads)) + ' reads loaded\n\n', flush=True, file=print_dest)
    return reads, check_reads, read_type


def get_albacore_barcode_from_path(albacore_path):
    if '/unclassified/' in albacore_path:
        return 'none'
    matches = re.findall('/barcode(\\d\\d)/', albacore_path)
    if matches:
        albacore_barcode_num = matches[-1]
        return 'BC' + albacore_barcode_num
    return None


def find_matching_adapter_sets(check_reads, verbosity, end_size, scoring_scheme_vals, print_dest,
                               adapter_threshold, threads):
    """
    Aligns all of the adapter sets to the start/end of reads to see which (if any) matches best.
    """
    read_count = len(check_reads)
    if verbosity > 0:
        print(bold_underline('Looking for known adapter sets'), flush=True, file=print_dest)
        output_progress_line(0, read_count, print_dest)

    search_adapters = [a for a in ADAPTERS if '(full sequence)' not in a.name]
    search_adapter_count = len(search_adapters)

    # If single-threaded, do the work in a simple loop.
    if threads == 1:
        for read_num, read in enumerate(check_reads):
            for adapter_set in search_adapters:
                read.align_adapter_set(adapter_set, end_size, scoring_scheme_vals)
            if verbosity > 0:
                output_progress_line(read_num+1, read_count, print_dest)

    # If multi-threaded, use a thread pool.
    else:
        def align_adapter_set_one_arg(all_args):
            r, a, b, c = all_args
            r.align_adapter_set(a, b, c)
        with ThreadPool(threads) as pool:
            arg_list = []
            for read in check_reads:
                for adapter_set in search_adapters:
                    arg_list.append((read, adapter_set, end_size, scoring_scheme_vals))
            finished_count = 0
            for _ in pool.imap(align_adapter_set_one_arg, arg_list):
                finished_count += 1
                if verbosity > 0 and finished_count % search_adapter_count == 0:
                    output_progress_line(finished_count // search_adapter_count,
                                         read_count, print_dest)

    if verbosity > 0:
        output_progress_line(read_count, read_count, print_dest, end_newline=True)

    return [x for x in search_adapters if x.best_start_or_end_score() >= adapter_threshold]


def choose_barcoding_kit(adapter_sets, verbosity, print_dest):
    """
    If the user is sorting reads by barcode bin, choose one barcode configuration (rev comp
    barcodes at the start of the read or at the end of the read) and ignore the other.
    """
    # Tally up scores for forward and reverse barcodes.
    forward_start_or_end, reverse_start_or_end = 0, 0
    forward_start_and_end, reverse_start_and_end = 0, 0
    for adapter_set in adapter_sets:
        if 'barcode' in adapter_set.name.lower():
            if '(forward)' in adapter_set.name.lower():
                forward_start_or_end += adapter_set.best_start_or_end_score()
                forward_start_and_end += adapter_set.best_start_score
                forward_start_and_end += adapter_set.best_end_score
            elif '(reverse)' in adapter_set.name.lower():
                reverse_start_or_end += adapter_set.best_start_or_end_score()
                reverse_start_and_end += adapter_set.best_start_score
                reverse_start_and_end += adapter_set.best_end_score

    if forward_start_or_end == 0 and reverse_start_or_end == 0:
        sys.exit('Error: no barcodes were found, so Porechop cannot perform barcode demultiplexing')

    # If possible, make a decision using each barcode's best start OR end score.
    orientation = None
    if forward_start_or_end > reverse_start_or_end:
        orientation = 'forward'
    elif reverse_start_or_end > forward_start_or_end:
        orientation = 'reverse'

    # If that didn't work (i.e. it's a tie between forward and reverse), then choose based on the
    # sum of both start AND end scores.
    elif forward_start_and_end > reverse_start_and_end:
        orientation = 'forward'
    elif reverse_start_and_end > forward_start_and_end:
        orientation = 'reverse'

    if orientation is None:
        sys.exit('Error: Porechop could not determine barcode orientation')

    if verbosity > 0:
        print('\nBarcodes determined to be in ' + orientation + ' orientation', file=print_dest)
    return orientation


def fix_up_1d2_sets(matching_sets):
    """
    The 1D^2 adapters share a lot of common sequence with the old SQK-MAP006_short adapters, so if
    there are strong 1D^2 hits, we can exclude the redundant SQK-MAP006_short adapters.
    """
    matching_set_names = [x.name for x in matching_sets]
    if '1D^2 part 1' in matching_set_names and '1D^2 part 2' in matching_set_names and \
            'SQK-MAP006 Short' in matching_set_names:
        part_1_score = [x for x in matching_sets
                        if x.name == '1D^2 part 1'][0].best_start_or_end_score()
        part_2_score = [x for x in matching_sets
                        if x.name == '1D^2 part 2'][0].best_start_or_end_score()
        sqk_score = [x for x in matching_sets
                     if x.name == 'SQK-MAP006 Short'][0].best_start_or_end_score()
        if part_1_score >= sqk_score and part_2_score >= sqk_score:
            matching_sets = [x for x in matching_sets if x.name != 'SQK-MAP006 Short']
    return matching_sets


def display_adapter_set_results(matching_sets, verbosity, print_dest):
    if verbosity < 1:
        return
    table = [['Set', 'Best read start %ID', 'Best read end %ID']]
    row_colours = {}
    matching_set_names = [x.name for x in matching_sets]
    search_adapters = [a for a in ADAPTERS if '(full sequence)' not in a.name]
    for adapter_set in search_adapters:
        start_score = '%.1f' % adapter_set.best_start_score
        end_score = '%.1f' % adapter_set.best_end_score
        table.append([adapter_set.name, start_score, end_score])
        if adapter_set.name in matching_set_names:
            row_colours[len(table) - 1] = 'green'
    if verbosity > 0:
        print_table(table, print_dest, alignments='LRR', row_colour=row_colours,
                    fixed_col_widths=[35, 8, 8])


def add_full_barcode_adapter_sets(matching_sets):
    """
    This function adds some new 'full' adapter sequences based on what was already found. For
    example, if the ligation adapters and the reverse barcode adapters are found, it assumes we are
    looking at a native barcoding run and so it adds the complete native barcoding adapter
    sequences (with the barcode's upstream and downstream context included).
    """
    matching_set_names = [x.name for x in matching_sets]

    for i in range(1, 97):

        # Native barcode full sequences
        if all(x in matching_set_names
               for x in ['SQK-NSK007', 'Barcode ' + str(i) + ' (reverse)']):
            matching_sets.append(make_full_native_barcode_adapter(i))

        # Rapid barcode full sequences
        if all(x in matching_set_names
               for x in ['Rapid', 'Barcode ' + str(i) + ' (forward)']):
            if 'RBK004_upstream' in matching_set_names:
                matching_sets.append(make_new_full_rapid_barcode_adapter(i))
            elif 'SQK-NSK007' in matching_set_names:
                matching_sets.append(make_old_full_rapid_barcode_adapter(i))

    return matching_sets


def find_adapters_at_read_ends(reads, matching_sets, verbosity, end_size, extra_trim_size,
                               end_threshold, scoring_scheme_vals, print_dest, min_trim_size,
                               threads, check_barcodes, barcode_threshold, barcode_diff,
                               require_two_barcodes, forward_or_reverse_barcodes):
    if verbosity > 0:
        print(bold_underline('Trimming adapters from read ends'),
              file=print_dest)
        name_len = max(max(len(x.start_sequence[0])
                           if x.start_sequence else 0 for x in matching_sets),
                       max(len(x.end_sequence[0])
                           if x.end_sequence else 0 for x in matching_sets))
        for matching_set in matching_sets:
            if matching_set.start_sequence:
                print('  ' + matching_set.start_sequence[0].rjust(name_len) + ': ' +
                      red(matching_set.start_sequence[1]), file=print_dest)
            if matching_set.end_sequence:
                print('  ' + matching_set.end_sequence[0].rjust(name_len) + ': ' +
                      red(matching_set.end_sequence[1]), file=print_dest)
        print('', file=print_dest)

    read_count = len(reads)
    if verbosity == 1:
        output_progress_line(0, read_count, print_dest)

    # If single-threaded, do the work in a simple loop.
    if threads == 1:
        for read_num, read in enumerate(reads):
            read.find_start_trim(matching_sets, end_size, extra_trim_size, end_threshold,
                                 scoring_scheme_vals, min_trim_size, check_barcodes,
                                 forward_or_reverse_barcodes)
            read.find_end_trim(matching_sets, end_size, extra_trim_size, end_threshold,
                               scoring_scheme_vals, min_trim_size, check_barcodes,
                               forward_or_reverse_barcodes)
            if check_barcodes:
                read.determine_barcode(barcode_threshold, barcode_diff, require_two_barcodes)
            if verbosity == 1:
                output_progress_line(read_num+1, read_count, print_dest)
            elif verbosity == 2:
                print(read.formatted_start_and_end_seq(end_size, extra_trim_size, check_barcodes),
                      file=print_dest)
            elif verbosity > 2:
                print(read.full_start_end_output(end_size, extra_trim_size, check_barcodes),
                      file=print_dest)

    # If multi-threaded, use a thread pool.
    else:
        def start_end_trim_one_arg(all_args):
            r, a, b, c, d, e, f, g, h, i, j, k, v = all_args
            r.find_start_trim(a, b, c, d, e, f, g, k)
            r.find_end_trim(a, b, c, d, e, f, g, k)
            if check_barcodes:
                r.determine_barcode(h, i, j)
            if v == 2:
                return r.formatted_start_and_end_seq(b, c, g)
            if v > 2:
                return r.full_start_end_output(b, c, g)
            else:
                return ''
        with ThreadPool(threads) as pool:
            arg_list = []
            for read in reads:
                arg_list.append((read, matching_sets, end_size, extra_trim_size, end_threshold,
                                 scoring_scheme_vals, min_trim_size, check_barcodes,
                                 barcode_threshold, barcode_diff, require_two_barcodes,
                                 forward_or_reverse_barcodes, verbosity))
            finished_count = 0
            for out in pool.imap(start_end_trim_one_arg, arg_list):
                finished_count += 1
                if verbosity == 1:
                    output_progress_line(finished_count, read_count, print_dest)
                elif verbosity > 1:
                    print(out, file=print_dest, flush=True)

    if verbosity == 1:
        output_progress_line(read_count, read_count, print_dest, end_newline=True)
    if verbosity > 0:
        print('', file=print_dest)


def display_read_end_trimming_summary(reads, verbosity, print_dest):
    if verbosity < 1:
        return
    start_trim_total = sum(x.start_trim_amount for x in reads)
    start_trim_count = sum(1 if x.start_trim_amount else 0 for x in reads)
    end_trim_count = sum(1 if x.end_trim_amount else 0 for x in reads)
    end_trim_total = sum(x.end_trim_amount for x in reads)
    print(int_to_str(start_trim_count).rjust(len(int_to_str(len(reads)))) + ' / ' +
          int_to_str(len(reads)) + ' reads had adapters trimmed from their start (' +
          int_to_str(start_trim_total) + ' bp removed)', file=print_dest)
    print(int_to_str(end_trim_count).rjust(len(int_to_str(len(reads)))) + ' / ' +
          int_to_str(len(reads)) + ' reads had adapters trimmed from their end (' +
          int_to_str(end_trim_total) + ' bp removed)', file=print_dest)
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
        if matching_set.start_sequence:
            adapters.append(matching_set.start_sequence)
        if matching_set.end_sequence:
            if (not matching_set.start_sequence) or \
                    matching_set.end_sequence[1] != matching_set.start_sequence[1]:
                adapters.append(matching_set.end_sequence)

    start_sequence_names = set()
    end_sequence_names = set()
    for matching_set in matching_sets:
        if matching_set.start_sequence:
            start_sequence_names.add(matching_set.start_sequence[0])
        if matching_set.end_sequence:
            end_sequence_names.add(matching_set.end_sequence[0])

    read_count = len(reads)
    if verbosity == 1:
        output_progress_line(0, read_count, print_dest)

    # If single-threaded, do the work in a simple loop.
    if threads == 1:
        for read_num, read in enumerate(reads):
            read.find_middle_adapters(adapters, middle_threshold, extra_trim_good_side,
                                      extra_trim_bad_side, scoring_scheme_vals,
                                      start_sequence_names, end_sequence_names)
            if verbosity == 1:
                output_progress_line(read_num+1, read_count, print_dest)
            if read.middle_adapter_positions and verbosity > 1:
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
            finished_count = 0
            for out in pool.imap(find_middle_adapters_one_arg, arg_list):
                finished_count += 1
                if verbosity == 1:
                    output_progress_line(finished_count + 1, read_count, print_dest)
                if verbosity > 1 and out:
                    print(out, file=print_dest, flush=True)

    if verbosity == 1:
        output_progress_line(read_count, read_count, print_dest, end_newline=True)
        print('', flush=True, file=print_dest)


def display_read_middle_trimming_summary(reads, discard_middle, verbosity, print_dest):
    if verbosity < 1:
        return
    middle_trim_count = sum(1 if x.middle_adapter_positions else 0 for x in reads)
    verb = 'discarded' if discard_middle else 'split'
    print(int_to_str(middle_trim_count) + ' / ' + int_to_str(len(reads)) + ' reads were ' + verb +
          ' based on middle adapters\n\n', file=print_dest)


def output_reads(reads, out_format, output, read_type, verbosity, discard_middle,
                 min_split_size, print_dest, barcode_dir, input_filename,
                 untrimmed, threads, discard_unassigned):
    if verbosity > 0:
        trimmed_or_untrimmed = 'untrimmed' if untrimmed else 'trimmed'
        if barcode_dir is not None:
            verb = 'Saving '
            destination = 'barcode-specific files'
        elif output is None:
            verb = 'Outputting '
            destination = 'stdout'
        else:
            verb = 'Saving '
            destination = 'file'
        print(bold_underline(verb + trimmed_or_untrimmed + ' reads to ' + destination),
              flush=True, file=print_dest)

    if out_format == 'auto':
        if output is None:
            out_format = read_type.lower()
            if barcode_dir is not None and input_filename.lower().endswith('.gz'):
                out_format += '.gz'
        elif '.fasta.gz' in output.lower():
            out_format = 'fasta.gz'
        elif '.fastq.gz' in output.lower():
            out_format = 'fastq.gz'
        elif '.fasta' in output.lower():
            out_format = 'fasta'
        elif '.fastq' in output.lower():
            out_format = 'fastq'
        else:
            out_format = read_type.lower()

    gzipped_out = False
    gzip_command = 'gzip'
    if out_format.endswith('.gz') and (barcode_dir is not None or output is not None):
        gzipped_out = True
        out_format = out_format[:-3]
        if shutil.which('pigz'):
            if verbosity > 0:
                print('pigz found - using it to compress instead of gzip')
            gzip_command = 'pigz -p ' + str(threads)
        else:
            if verbosity > 0:
                print('pigz not found - using gzip to compress')

    # Output reads to barcode bins.
    if barcode_dir is not None:
        if not os.path.isdir(barcode_dir):
            os.makedirs(barcode_dir)
        barcode_files = {}
        barcode_read_counts, barcode_base_counts = defaultdict(int), defaultdict(int)
        for read in reads:
            barcode_name = read.barcode_call
            if discard_unassigned and barcode_name == 'none':
                continue
            if out_format == 'fasta':
                read_str = read.get_fasta(min_split_size, discard_middle, untrimmed)
            else:
                read_str = read.get_fastq(min_split_size, discard_middle, untrimmed)
            if not read_str:
                continue
            if barcode_name not in barcode_files:
                barcode_files[barcode_name] = \
                    open(os.path.join(barcode_dir, barcode_name + '.' + out_format), 'wt')
            barcode_files[barcode_name].write(read_str)
            barcode_read_counts[barcode_name] += 1
            if untrimmed:
                seq_length = len(read.seq)
            else:
                seq_length = read.seq_length_with_start_end_adapters_trimmed()
            barcode_base_counts[barcode_name] += seq_length
        table = [['Barcode', 'Reads', 'Bases', 'File']]

        for barcode_name in sorted(barcode_files.keys()):
            barcode_files[barcode_name].close()
            bin_filename = os.path.join(barcode_dir, barcode_name + '.' + out_format)

            if gzipped_out:
                if not os.path.isfile(bin_filename):
                    continue
                bin_filename_gz = bin_filename + '.gz'
                if os.path.isfile(bin_filename_gz):
                    os.remove(bin_filename_gz)
                try:
                    subprocess.check_output(gzip_command + ' ' + bin_filename,
                                            stderr=subprocess.STDOUT, shell=True)
                except subprocess.CalledProcessError:
                    pass
                bin_filename = bin_filename_gz

            table_row = [barcode_name, int_to_str(barcode_read_counts[barcode_name]),
                         int_to_str(barcode_base_counts[barcode_name]), bin_filename]
            table.append(table_row)

        if verbosity > 0:
            print('')
            print_table(table, print_dest, alignments='LRRL', max_col_width=60, col_separation=2)

    # Output to all reads to stdout.
    elif output is None:
        for read in reads:
            read_str = read.get_fasta(min_split_size, discard_middle) if out_format == 'fasta' \
                else read.get_fastq(min_split_size, discard_middle)
            print(read_str, end='')
        if verbosity > 0:
            print('Done', flush=True, file=print_dest)

    # Output to all reads to file.
    else:
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
            subprocess.check_output(gzip_command + ' -c ' + out_filename + ' > ' + output,
                                    stderr=subprocess.STDOUT, shell=True)
            os.remove(out_filename)
        if verbosity > 0:
            print('\nSaved result to ' + os.path.abspath(output), file=print_dest)

    if verbosity > 0:
        print('', flush=True, file=print_dest)


def output_progress_line(completed, total, print_dest, end_newline=False, step=10):
    if step > 1 and completed % step != 0 and completed != total:
        return
    progress_str = int_to_str(completed) + ' / ' + int_to_str(total)
    if total > 0:
        percent = 100.0 * completed / total
    else:
        percent = 0.0
    progress_str += ' (' + '%.1f' % percent + '%)'

    end_char = '\n' if end_newline else ''
    print('\r' + progress_str, end=end_char, flush=True, file=print_dest)
