"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Porechop

This module contains some tests for Porechop. To run them, execute `python3 -m unittest` from the
root Porechop directory.

This file is part of Porechop. Porechop is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Porechop is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Porechop. If
not, see <http://www.gnu.org/licenses/>.
"""

import unittest
import os
import subprocess
import shutil
import porechop.misc


def get_read_type(filename):
    _, read_type = porechop.misc.load_fasta_or_fastq(filename)
    read_type = read_type.lower()
    compression_type = porechop.misc.get_compression_type(filename)
    if compression_type == 'plain':
        return read_type
    elif compression_type == 'gz':
        return read_type + '.gz'
    else:
        assert False


class TestOutputFormat(unittest.TestCase):

    def run_command(self, command, input_filename):
        runner_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'porechop-runner.py')
        input_path = os.path.join(os.path.dirname(__file__), input_filename)
        command = command.replace('porechop', runner_path)
        command = command.replace('IN', input_path)
        output_name = 'TEMP_' + str(os.getpid())
        command = command.replace('OUT', output_name)
        try:
            self.output_file = [x for x in command.split() if output_name in x][0]
        except IndexError:
            self.output_file = ''
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = p.communicate()
        return out.decode(), err.decode()

    def tearDown(self):
        if os.path.isfile(self.output_file):
            os.remove(self.output_file)

    # The following tests use the auto format and should determine the output type using the output
    # filename.

    def test_auto_format_fastq_to_fastq(self):
        self.run_command('porechop -i IN -o OUT.fastq', 'test_format.fastq')
        self.assertEqual(get_read_type(self.output_file), 'fastq')

    def test_auto_format_fastq_to_fastq_gz(self):
        self.run_command('porechop -i IN -o OUT.fastq.gz', 'test_format.fastq')
        self.assertEqual(get_read_type(self.output_file), 'fastq.gz')

    def test_auto_format_fastq_to_fasta(self):
        self.run_command('porechop -i IN -o OUT.fasta', 'test_format.fastq')
        self.assertEqual(get_read_type(self.output_file), 'fasta')

    def test_auto_format_fastq_to_fasta_gz(self):
        self.run_command('porechop -i IN -o OUT.fasta.gz', 'test_format.fastq')
        self.assertEqual(get_read_type(self.output_file), 'fasta.gz')

    def test_auto_format_fastq_to_fastq_gz_weird_caps(self):
        self.run_command('porechop -i IN -o OUT.fAsTq.Gz', 'test_format.fastq')
        self.assertEqual(get_read_type(self.output_file), 'fastq.gz')

    # The following tests use the auto format and pipe to output. They determine the output format
    # from the input filename, but won't gzip the reads.

    def test_auto_format_fastq_to_pipe(self):
        self.run_command('porechop -i IN > OUT', 'test_format.fastq')
        self.assertEqual(get_read_type(self.output_file), 'fastq')

    def test_auto_format_fastq_gz_to_pipe(self):
        self.run_command('porechop -i IN > OUT', 'test_format.fastq.gz')
        self.assertEqual(get_read_type(self.output_file), 'fastq')

    def test_auto_format_fasta_to_pipe(self):
        self.run_command('porechop -i IN > OUT', 'test_format.fasta')
        self.assertEqual(get_read_type(self.output_file), 'fasta')

    def test_auto_format_fasta_gz_to_pipe(self):
        self.run_command('porechop -i IN > OUT', 'test_format.fasta.gz')
        self.assertEqual(get_read_type(self.output_file), 'fasta')

    # The following tests explicitly set the output format.

    def test_explicit_format_fastq_to_fastq(self):
        self.run_command('porechop -i IN -o OUT --format fastq', 'test_format.fastq')
        self.assertEqual(get_read_type(self.output_file), 'fastq')

    def test_explicit_format_fastq_to_fastq_gz(self):
        self.run_command('porechop -i IN -o OUT --format fastq.gz', 'test_format.fastq')
        self.assertEqual(get_read_type(self.output_file), 'fastq.gz')

    def test_explicit_format_fastq_to_fasta(self):
        self.run_command('porechop -i IN -o OUT --format fasta', 'test_format.fastq')
        self.assertEqual(get_read_type(self.output_file), 'fasta')

    def test_explicit_format_fastq_to_fasta_gz(self):
        self.run_command('porechop -i IN -o OUT --format fasta.gz', 'test_format.fastq')
        self.assertEqual(get_read_type(self.output_file), 'fasta.gz')

    def test_explicit_format_fasta_to_fasta(self):
        self.run_command('porechop -i IN -o OUT --format fasta', 'test_format.fasta')
        self.assertEqual(get_read_type(self.output_file), 'fasta')

    def test_explicit_format_fasta_to_fasta_gz(self):
        self.run_command('porechop -i IN -o OUT --format fasta.gz', 'test_format.fasta')
        self.assertEqual(get_read_type(self.output_file), 'fasta.gz')

    # For these tests the output filename conflicts with the requested format. The requested format
    # should be followed, not the filename.

    def test_explicit_format_fastq_to_fastq_conflicting_name_1(self):
        self.run_command('porechop -i IN -o OUT.fasta --format fastq', 'test_format.fastq')
        self.assertEqual(get_read_type(self.output_file), 'fastq')

    def test_explicit_format_fastq_to_fastq_conflicting_name_2(self):
        self.run_command('porechop -i IN -o OUT.fasta.gz --format fastq', 'test_format.fastq.gz')
        self.assertEqual(get_read_type(self.output_file), 'fastq')

    def test_explicit_format_fastq_to_fasta_gz_conflicting_name_1(self):
        self.run_command('porechop -i IN -o OUT.fastq --format fasta.gz', 'test_format.fastq')
        self.assertEqual(get_read_type(self.output_file), 'fasta.gz')

    def test_explicit_format_fastq_to_fasta_gz_conflicting_name_2(self):
        self.run_command('porechop -i IN -o OUT.fastq.gz --format fasta.gz', 'test_format.fasta')
        self.assertEqual(get_read_type(self.output_file), 'fasta.gz')


class TestOutputFormatBarcodes(unittest.TestCase):

    def run_command(self, command, input_filename):
        runner_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'porechop-runner.py')
        input_path = os.path.join(os.path.dirname(__file__), input_filename)
        command = command.replace('porechop', runner_path)
        command = command.replace('IN', input_path)
        self.output_dir = 'TEMP_' + str(os.getpid())
        command = command.replace('B_DIR', self.output_dir)
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = p.communicate()
        return out.decode(), err.decode()

    def tearDown(self):
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)

    def check_output_file_formats(self, output_files, expected_format):
        for output_file in [os.path.join(self.output_dir, x) for x in output_files]:
            self.assertTrue(os.path.isfile(output_file))
            self.assertEqual(get_read_type(output_file), expected_format)

    # The following tests use the auto format and should determine the output type using the input
    # type.
    
    def test_auto_format_fastq_to_fastq(self):
        self.run_command('porechop -i IN -b B_DIR', 'test_format_barcodes.fastq')
        self.check_output_file_formats(['BC01.fastq', 'BC02.fastq', 'BC03.fastq', 'none.fastq'],
                                       'fastq')

    def test_auto_format_fastq_gz_to_fastq_gz(self):
        self.run_command('porechop -i IN -b B_DIR', 'test_format_barcodes.fastq.gz')
        self.check_output_file_formats(['BC01.fastq.gz', 'BC02.fastq.gz', 'BC03.fastq.gz',
                                        'none.fastq.gz'], 'fastq.gz')

    def test_auto_format_fasta_to_fasta(self):
        self.run_command('porechop -i IN -b B_DIR', 'test_format_barcodes.fasta')
        self.check_output_file_formats(['BC01.fasta', 'BC02.fasta', 'BC03.fasta', 'none.fasta'],
                                       'fasta')

    def test_auto_format_fasta_gz_to_fasta_gz(self):
        self.run_command('porechop -i IN -b B_DIR', 'test_format_barcodes.fasta.gz')
        self.check_output_file_formats(['BC01.fasta.gz', 'BC02.fasta.gz', 'BC03.fasta.gz',
                                        'none.fasta.gz'], 'fasta.gz')

    # The following tests explicitly set the output format.

    def test_explicit_format_fastq_to_fastq(self):
        self.run_command('porechop -i IN -b B_DIR --format fastq', 'test_format_barcodes.fastq')
        self.check_output_file_formats(['BC01.fastq', 'BC02.fastq', 'BC03.fastq', 'none.fastq'],
                                       'fastq')

    def test_explicit_format_fastq_to_fastq_gz(self):
        self.run_command('porechop -i IN -b B_DIR --format fastq.gz', 'test_format_barcodes.fastq')
        self.check_output_file_formats(['BC01.fastq.gz', 'BC02.fastq.gz', 'BC03.fastq.gz',
                                        'none.fastq.gz'], 'fastq.gz')

    def test_explicit_format_fastq_gz_to_fastq(self):
        self.run_command('porechop -i IN -b B_DIR --format fastq', 'test_format_barcodes.fastq.gz')
        self.check_output_file_formats(['BC01.fastq', 'BC02.fastq', 'BC03.fastq', 'none.fastq'],
                                       'fastq')

    def test_explicit_format_fastq_gz_to_fastq_gz(self):
        self.run_command('porechop -i IN -b B_DIR --format fastq.gz',
                         'test_format_barcodes.fastq.gz')
        self.check_output_file_formats(['BC01.fastq.gz', 'BC02.fastq.gz', 'BC03.fastq.gz',
                                        'none.fastq.gz'], 'fastq.gz')

    def test_explicit_format_fastq_to_fasta(self):
        self.run_command('porechop -i IN -b B_DIR --format fasta', 'test_format_barcodes.fastq')
        self.check_output_file_formats(['BC01.fasta', 'BC02.fasta', 'BC03.fasta', 'none.fasta'],
                                       'fasta')

    def test_explicit_format_fastq_to_fasta_gz(self):
        self.run_command('porechop -i IN -b B_DIR --format fasta.gz', 'test_format_barcodes.fastq')
        self.check_output_file_formats(['BC01.fasta.gz', 'BC02.fasta.gz', 'BC03.fasta.gz',
                                        'none.fasta.gz'], 'fasta.gz')

    def test_explicit_format_fasta_to_fasta(self):
        self.run_command('porechop -i IN -b B_DIR --format fasta', 'test_format_barcodes.fasta')
        self.check_output_file_formats(['BC01.fasta', 'BC02.fasta', 'BC03.fasta', 'none.fasta'],
                                       'fasta')

    def test_explicit_format_fasta_to_fasta_gz(self):
        self.run_command('porechop -i IN -b B_DIR --format fasta.gz', 'test_format_barcodes.fasta')
        self.check_output_file_formats(['BC01.fasta.gz', 'BC02.fasta.gz', 'BC03.fasta.gz',
                                        'none.fasta.gz'], 'fasta.gz')

    def test_explicit_format_fasta_gz_to_fasta(self):
        self.run_command('porechop -i IN -b B_DIR --format fasta', 'test_format_barcodes.fasta.gz')
        self.check_output_file_formats(['BC01.fasta', 'BC02.fasta', 'BC03.fasta', 'none.fasta'],
                                       'fasta')

    def test_explicit_format_fasta_gz_to_fasta_gz(self):
        self.run_command('porechop -i IN -b B_DIR --format fasta.gz',
                         'test_format_barcodes.fasta.gz')
        self.check_output_file_formats(['BC01.fasta.gz', 'BC02.fasta.gz', 'BC03.fasta.gz',
                                        'none.fasta.gz'], 'fasta.gz')
