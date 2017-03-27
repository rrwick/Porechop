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
import shutil
import subprocess
import porechop.misc


class TestBarcodes(unittest.TestCase):
    """
    Tests barcode demultiplexing.
    """
    def run_command(self, command):
        runner_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'porechop-runner.py')
        input_path = os.path.join(os.path.dirname(__file__), 'test_barcodes.fastq')
        command = command.replace('porechop', runner_path)
        command = command.replace('INPUT', input_path)
        self.output_dir = 'TEMP_' + str(os.getpid())
        command = command.replace('BARCODE_DIR', self.output_dir)
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = p.communicate()
        return out.decode(), err.decode()

    def tearDown(self):
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)

    def load_trimmed_reads(self, filename):
        full_filepath = os.path.join(self.output_dir, filename)
        trimmed_reads, read_type = porechop.misc.load_fasta_or_fastq(full_filepath)
        if read_type == 'FASTA':
            trimmed_reads = [(x[2], x[1], '') for x in trimmed_reads]
        else:  # FASTQ
            trimmed_reads = [(x[4], x[1], x[3]) for x in trimmed_reads]
        return trimmed_reads

    def test_barcodes_1(self):
        """
        Tests with default settings.
        """
        out, _ = self.run_command('porechop -i INPUT -b BARCODE_DIR')

        nb01_trimmed_reads = self.load_trimmed_reads('NB01.fastq.gz')
        nb02_trimmed_reads = self.load_trimmed_reads('NB02.fastq.gz')
        nb03_trimmed_reads = self.load_trimmed_reads('NB03.fastq.gz')
        none_trimmed_reads = self.load_trimmed_reads('none.fastq.gz')

        nb01_read_names = sorted(x[0] for x in nb01_trimmed_reads)
        nb02_read_names = sorted(x[0] for x in nb02_trimmed_reads)
        nb03_read_names = sorted(x[0] for x in nb03_trimmed_reads)
        none_read_names = sorted(x[0] for x in none_trimmed_reads)

        self.assertEqual(nb01_read_names, ['1', '4'])
        self.assertEqual(nb02_read_names, ['2', '5'])
        self.assertEqual(nb03_read_names, ['3'])
        self.assertEqual(none_read_names, ['6', '8'])

        self.assertEqual(sum(len(x[1]) for x in nb01_trimmed_reads), 8994)
        self.assertEqual(sum(len(x[2]) for x in nb01_trimmed_reads), 8994)
        self.assertTrue('NB01         2   8,994' in out)

        self.assertEqual(sum(len(x[1]) for x in nb02_trimmed_reads), 9394)
        self.assertEqual(sum(len(x[2]) for x in nb02_trimmed_reads), 9394)
        self.assertTrue('NB02         2   9,394' in out)

        self.assertEqual(sum(len(x[1]) for x in nb03_trimmed_reads), 6996)
        self.assertEqual(sum(len(x[2]) for x in nb03_trimmed_reads), 6996)
        self.assertTrue('NB03         1   6,996' in out)

        self.assertEqual(sum(len(x[1]) for x in none_trimmed_reads), 13496)
        self.assertEqual(sum(len(x[2]) for x in none_trimmed_reads), 13496)
        self.assertTrue('none         2  13,496' in out)

    def test_barcodes_2(self):
        """
        Tests with --require_two_barcodes.
        """
        out, _ = self.run_command('porechop -i INPUT -b BARCODE_DIR --require_two_barcodes')

        nb01_trimmed_reads = self.load_trimmed_reads('NB01.fastq.gz')
        nb02_trimmed_reads = self.load_trimmed_reads('NB02.fastq.gz')
        nb03_trimmed_reads = self.load_trimmed_reads('NB03.fastq.gz')
        none_trimmed_reads = self.load_trimmed_reads('none.fastq.gz')

        nb01_read_names = sorted(x[0] for x in nb01_trimmed_reads)
        nb02_read_names = sorted(x[0] for x in nb02_trimmed_reads)
        nb03_read_names = sorted(x[0] for x in nb03_trimmed_reads)
        none_read_names = sorted(x[0] for x in none_trimmed_reads)

        self.assertEqual(nb01_read_names, ['1'])
        self.assertEqual(nb02_read_names, ['2'])
        self.assertEqual(nb03_read_names, ['3'])
        self.assertEqual(none_read_names, ['4', '5', '6', '8'])

        self.assertEqual(sum(len(x[1]) for x in nb01_trimmed_reads), 4996)
        self.assertEqual(sum(len(x[2]) for x in nb01_trimmed_reads), 4996)
        self.assertTrue('NB01         1   4,996' in out)

        self.assertEqual(sum(len(x[1]) for x in nb02_trimmed_reads), 6096)
        self.assertEqual(sum(len(x[2]) for x in nb02_trimmed_reads), 6096)
        self.assertTrue('NB02         1   6,096' in out)

        self.assertEqual(sum(len(x[1]) for x in nb03_trimmed_reads), 6996)
        self.assertEqual(sum(len(x[2]) for x in nb03_trimmed_reads), 6996)
        self.assertTrue('NB03         1   6,996' in out)

        self.assertEqual(sum(len(x[1]) for x in none_trimmed_reads), 20792)
        self.assertEqual(sum(len(x[2]) for x in none_trimmed_reads), 20792)
        self.assertTrue('none         4  20,792' in out)
