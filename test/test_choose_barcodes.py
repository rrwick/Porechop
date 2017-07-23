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
    Tests barcode demultiplexing where there are conflicting forward and reverse barcodes. In such
    a case, Porechop should choose the more popular option.
    """
    def run_command(self, command, input_path):
        runner_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'porechop-runner.py')
        input_path = os.path.join(os.path.dirname(__file__), input_path)
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

    def test_choose_forward(self):
        """
        In this test, forward barcodes are more common and reverse barcodes should be ignored.
        """
        out, _ = self.run_command('porechop -i INPUT -b BARCODE_DIR',
                                  'test_choose_barcodes_1.fasta')

        bc01_trimmed_reads = self.load_trimmed_reads('BC01.fasta')
        bc02_trimmed_reads = self.load_trimmed_reads('BC02.fasta')
        bc03_trimmed_reads = self.load_trimmed_reads('BC03.fasta')
        bc04_trimmed_reads = self.load_trimmed_reads('BC04.fasta')
        bc05_trimmed_reads = self.load_trimmed_reads('BC05.fasta')
        bc06_trimmed_reads = self.load_trimmed_reads('BC06.fasta')
        none_trimmed_reads = self.load_trimmed_reads('none.fasta')

        bc01_read_names = sorted(x[0] for x in bc01_trimmed_reads)
        bc02_read_names = sorted(x[0] for x in bc02_trimmed_reads)
        bc03_read_names = sorted(x[0] for x in bc03_trimmed_reads)
        bc04_read_names = sorted(x[0] for x in bc04_trimmed_reads)
        bc05_read_names = sorted(x[0] for x in bc05_trimmed_reads)
        bc06_read_names = sorted(x[0] for x in bc06_trimmed_reads)
        none_read_names = sorted(x[0] for x in none_trimmed_reads)

        self.assertEqual(bc01_read_names, ['1'])
        self.assertEqual(bc02_read_names, ['2'])
        self.assertEqual(bc03_read_names, ['3'])
        self.assertEqual(bc04_read_names, ['4'])
        self.assertEqual(bc05_read_names, ['5'])
        self.assertEqual(bc06_read_names, ['6'])
        self.assertEqual(none_read_names, ['7', '8'])

        self.assertTrue('Barcodes determined to be in forward orientation' in out)
        self.assertTrue('0 / 8 reads were discarded based on middle adapters' in out)

    def test_choose_reverse(self):
        """
        In this test, reverse barcodes are more common and forward barcodes should be ignored.
        """
        out, _ = self.run_command('porechop -i INPUT -b BARCODE_DIR',
                                  'test_choose_barcodes_2.fasta')

        bc01_trimmed_reads = self.load_trimmed_reads('BC01.fasta')
        bc02_trimmed_reads = self.load_trimmed_reads('BC02.fasta')
        bc03_trimmed_reads = self.load_trimmed_reads('BC03.fasta')
        bc04_trimmed_reads = self.load_trimmed_reads('BC04.fasta')
        bc05_trimmed_reads = self.load_trimmed_reads('BC05.fasta')
        bc06_trimmed_reads = self.load_trimmed_reads('BC06.fasta')
        none_trimmed_reads = self.load_trimmed_reads('none.fasta')

        bc01_read_names = sorted(x[0] for x in bc01_trimmed_reads)
        bc02_read_names = sorted(x[0] for x in bc02_trimmed_reads)
        bc03_read_names = sorted(x[0] for x in bc03_trimmed_reads)
        bc04_read_names = sorted(x[0] for x in bc04_trimmed_reads)
        bc05_read_names = sorted(x[0] for x in bc05_trimmed_reads)
        bc06_read_names = sorted(x[0] for x in bc06_trimmed_reads)
        none_read_names = sorted(x[0] for x in none_trimmed_reads)

        self.assertEqual(bc01_read_names, ['1'])
        self.assertEqual(bc02_read_names, ['2'])
        self.assertEqual(bc03_read_names, ['3'])
        self.assertEqual(bc04_read_names, ['4'])
        self.assertEqual(bc05_read_names, ['5'])
        self.assertEqual(bc06_read_names, ['6'])
        self.assertEqual(none_read_names, ['7', '8'])

        self.assertTrue('Barcodes determined to be in reverse orientation' in out)
        self.assertTrue('0 / 8 reads were discarded based on middle adapters' in out)
