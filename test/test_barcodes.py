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

    def count_output_fastq_files(self):
        return len([x for x in os.listdir(self.output_dir)
                    if os.path.isfile(os.path.join(self.output_dir, x)) and 'fastq' in x])

    def test_barcodes_1(self):
        """
        Tests with default settings.
        """
        out, _ = self.run_command('porechop -i INPUT -b BARCODE_DIR')

        self.assertEqual(self.count_output_fastq_files(), 4)
        bc01_trimmed_reads = self.load_trimmed_reads('BC01.fastq')
        bc02_trimmed_reads = self.load_trimmed_reads('BC02.fastq')
        bc03_trimmed_reads = self.load_trimmed_reads('BC03.fastq')
        none_trimmed_reads = self.load_trimmed_reads('none.fastq')

        bc01_read_names = sorted(x[0] for x in bc01_trimmed_reads)
        bc02_read_names = sorted(x[0] for x in bc02_trimmed_reads)
        bc03_read_names = sorted(x[0] for x in bc03_trimmed_reads)
        none_read_names = sorted(x[0] for x in none_trimmed_reads)

        self.assertEqual(bc01_read_names, ['1', '4'])
        self.assertEqual(bc02_read_names, ['2', '5'])
        self.assertEqual(bc03_read_names, ['3'])
        self.assertEqual(none_read_names, ['6', '8'])

        self.assertEqual(sum(len(x[1]) for x in bc01_trimmed_reads), 8994)
        self.assertEqual(sum(len(x[2]) for x in bc01_trimmed_reads), 8994)
        self.assertTrue('BC01         2   8,994' in out)

        self.assertEqual(sum(len(x[1]) for x in bc02_trimmed_reads), 9394)
        self.assertEqual(sum(len(x[2]) for x in bc02_trimmed_reads), 9394)
        self.assertTrue('BC02         2   9,394' in out)

        self.assertEqual(sum(len(x[1]) for x in bc03_trimmed_reads), 6996)
        self.assertEqual(sum(len(x[2]) for x in bc03_trimmed_reads), 6996)
        self.assertTrue('BC03         1   6,996' in out)

        self.assertEqual(sum(len(x[1]) for x in none_trimmed_reads), 13496)
        self.assertEqual(sum(len(x[2]) for x in none_trimmed_reads), 13496)
        self.assertTrue('none         2  13,496' in out)

        self.assertTrue('Saving trimmed reads' in out)

    def test_barcodes_2(self):
        """
        Tests with --require_two_barcodes.
        """
        out, _ = self.run_command('porechop -i INPUT -b BARCODE_DIR --require_two_barcodes')

        self.assertEqual(self.count_output_fastq_files(), 4)
        bc01_trimmed_reads = self.load_trimmed_reads('BC01.fastq')
        bc02_trimmed_reads = self.load_trimmed_reads('BC02.fastq')
        bc03_trimmed_reads = self.load_trimmed_reads('BC03.fastq')
        none_trimmed_reads = self.load_trimmed_reads('none.fastq')

        bc01_read_names = sorted(x[0] for x in bc01_trimmed_reads)
        bc02_read_names = sorted(x[0] for x in bc02_trimmed_reads)
        bc03_read_names = sorted(x[0] for x in bc03_trimmed_reads)
        none_read_names = sorted(x[0] for x in none_trimmed_reads)

        self.assertEqual(bc01_read_names, ['1'])
        self.assertEqual(bc02_read_names, ['2'])
        self.assertEqual(bc03_read_names, ['3'])
        self.assertEqual(none_read_names, ['4', '5', '6', '8'])

        self.assertEqual(sum(len(x[1]) for x in bc01_trimmed_reads), 4996)
        self.assertEqual(sum(len(x[2]) for x in bc01_trimmed_reads), 4996)
        self.assertTrue('BC01         1   4,996' in out)

        self.assertEqual(sum(len(x[1]) for x in bc02_trimmed_reads), 6096)
        self.assertEqual(sum(len(x[2]) for x in bc02_trimmed_reads), 6096)
        self.assertTrue('BC02         1   6,096' in out)

        self.assertEqual(sum(len(x[1]) for x in bc03_trimmed_reads), 6996)
        self.assertEqual(sum(len(x[2]) for x in bc03_trimmed_reads), 6996)
        self.assertTrue('BC03         1   6,996' in out)

        self.assertEqual(sum(len(x[1]) for x in none_trimmed_reads), 20792)
        self.assertEqual(sum(len(x[2]) for x in none_trimmed_reads), 20792)
        self.assertTrue('none         4  20,792' in out)

        self.assertTrue('Saving trimmed reads' in out)

    def test_barcodes_3(self):
        """
        Tests with --untrimmed.
        """
        out, _ = self.run_command('porechop -i INPUT -b BARCODE_DIR --untrimmed')

        self.assertEqual(self.count_output_fastq_files(), 4)
        bc01_trimmed_reads = self.load_trimmed_reads('BC01.fastq')
        bc02_trimmed_reads = self.load_trimmed_reads('BC02.fastq')
        bc03_trimmed_reads = self.load_trimmed_reads('BC03.fastq')
        none_trimmed_reads = self.load_trimmed_reads('none.fastq')

        bc01_read_names = sorted(x[0] for x in bc01_trimmed_reads)
        bc02_read_names = sorted(x[0] for x in bc02_trimmed_reads)
        bc03_read_names = sorted(x[0] for x in bc03_trimmed_reads)
        none_read_names = sorted(x[0] for x in none_trimmed_reads)

        self.assertEqual(bc01_read_names, ['1', '4'])
        self.assertEqual(bc02_read_names, ['2', '5'])
        self.assertEqual(bc03_read_names, ['3'])
        self.assertEqual(none_read_names, ['6', '8'])

        self.assertEqual(sum(len(x[1]) for x in bc01_trimmed_reads), 9199)
        self.assertEqual(sum(len(x[2]) for x in bc01_trimmed_reads), 9199)
        self.assertTrue('BC01         2   9,199' in out)

        self.assertEqual(sum(len(x[1]) for x in bc02_trimmed_reads), 9594)
        self.assertEqual(sum(len(x[2]) for x in bc02_trimmed_reads), 9594)
        self.assertTrue('BC02         2   9,594' in out)

        self.assertEqual(sum(len(x[1]) for x in bc03_trimmed_reads), 7131)
        self.assertEqual(sum(len(x[2]) for x in bc03_trimmed_reads), 7131)
        self.assertTrue('BC03         1   7,131' in out)

        self.assertEqual(sum(len(x[1]) for x in none_trimmed_reads), 13631)
        self.assertEqual(sum(len(x[2]) for x in none_trimmed_reads), 13631)
        self.assertTrue('none         2  13,631' in out)

        self.assertTrue('Saving untrimmed reads' in out)

    def test_barcodes_4(self):
        """
        Tests with --discard_unassigned.
        """
        out, _ = self.run_command('porechop -i INPUT -b BARCODE_DIR --discard_unassigned')

        self.assertEqual(self.count_output_fastq_files(), 3)
        bc01_trimmed_reads = self.load_trimmed_reads('BC01.fastq')
        bc02_trimmed_reads = self.load_trimmed_reads('BC02.fastq')
        bc03_trimmed_reads = self.load_trimmed_reads('BC03.fastq')

        bc01_read_names = sorted(x[0] for x in bc01_trimmed_reads)
        bc02_read_names = sorted(x[0] for x in bc02_trimmed_reads)
        bc03_read_names = sorted(x[0] for x in bc03_trimmed_reads)

        self.assertEqual(bc01_read_names, ['1', '4'])
        self.assertEqual(bc02_read_names, ['2', '5'])
        self.assertEqual(bc03_read_names, ['3'])

        self.assertEqual(sum(len(x[1]) for x in bc01_trimmed_reads), 8994)
        self.assertEqual(sum(len(x[2]) for x in bc01_trimmed_reads), 8994)
        self.assertTrue('BC01         2  8,994' in out)

        self.assertEqual(sum(len(x[1]) for x in bc02_trimmed_reads), 9394)
        self.assertEqual(sum(len(x[2]) for x in bc02_trimmed_reads), 9394)
        self.assertTrue('BC02         2  9,394' in out)

        self.assertEqual(sum(len(x[1]) for x in bc03_trimmed_reads), 6996)
        self.assertEqual(sum(len(x[2]) for x in bc03_trimmed_reads), 6996)
        self.assertTrue('BC03         1  6,996' in out)

        self.assertTrue('Saving trimmed reads' in out)
