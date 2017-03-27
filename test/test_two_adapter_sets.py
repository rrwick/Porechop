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
import porechop.misc


class TestTwoAdapterSets(unittest.TestCase):
    """
    This test set contains a combination of two adapter sets: SQK-MAP006 and SQK-NSK007
    """
    def run_command(self, command):
        runner_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'porechop-runner.py')
        input_path = os.path.join(os.path.dirname(__file__), 'test_two_adapter_sets.fastq')
        command = command.replace('porechop', runner_path)
        command = command.replace('INPUT', input_path)
        output_name = 'TEMP_' + str(os.getpid())
        command = command.replace('OUTPUT', output_name)
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

    def load_trimmed_reads(self):
        trimmed_reads, read_type = porechop.misc.load_fasta_or_fastq(self.output_file)
        if read_type == 'FASTA':
            trimmed_reads = [(x[2], x[1], '') for x in trimmed_reads]
        else:  # FASTQ
            trimmed_reads = [(x[4], x[1], x[3]) for x in trimmed_reads]
        return trimmed_reads, read_type

    def check_trimmed_reads(self):
        """
        This function checks all of the trimmed reads (assuming default settings).
        """
        trimmed_reads, read_type = self.load_trimmed_reads()

        # Read 1 isn't trimmed.
        read_1 = trimmed_reads[0]
        self.assertEqual(read_1[0], '1')
        self.assertEqual(len(read_1[1]), 10000)
        self.assertTrue(read_1[1].startswith('GCGGAAGTCA'))
        self.assertTrue(read_1[1].endswith('TCCAGCCCCA'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_1[2]), 10000)
            self.assertTrue(read_1[2].startswith("4,)4+311.'"))
            self.assertTrue(read_1[2].endswith("3-16('...,"))

        # Read 2 has an adapter trimmed from the start.
        read_2 = trimmed_reads[1]
        self.assertEqual(read_2[0], '2')
        self.assertEqual(len(read_2[1]), 10970)
        self.assertTrue(read_2[1].startswith('TGCAGCACAA'))
        self.assertTrue(read_2[1].endswith('AAAGAACTGG'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_2[2]), 10970)
            self.assertTrue(read_2[2].startswith("&)*-++$+,."))
            self.assertTrue(read_2[2].endswith("(*'%)%.2+)"))

        # Read 3 has an adapter trimmed from the start.
        read_3 = trimmed_reads[2]
        self.assertEqual(read_3[0], '3')
        self.assertEqual(len(read_3[1]), 7970)
        self.assertTrue(read_3[1].startswith('TATTTCGCTT'))
        self.assertTrue(read_3[1].endswith('GAGAGTGTGA'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_3[2]), 7970)
            self.assertTrue(read_3[2].startswith("2:2.4%%%)$"))
            self.assertTrue(read_3[2].endswith("/-((+1-/%."))

        # Read 4 has two internal adapters (from different sets).
        read_4_1 = trimmed_reads[3]
        read_4_2 = trimmed_reads[4]
        read_4_3 = trimmed_reads[5]
        self.assertEqual(read_4_1[0], '4_1')
        self.assertEqual(read_4_2[0], '4_2')
        self.assertEqual(read_4_3[0], '4_3')
        self.assertEqual(len(read_4_1[1]), 1990)
        self.assertEqual(len(read_4_2[1]), 2868)
        self.assertEqual(len(read_4_3[1]), 1878)
        self.assertTrue(read_4_1[1].startswith('TCATGACAGT'))
        self.assertTrue(read_4_1[1].endswith('CGAAGTAGGG'))
        self.assertTrue(read_4_2[1].startswith('ACGTGCTCGA'))
        self.assertTrue(read_4_2[1].endswith('CCTTCCTAGA'))
        self.assertTrue(read_4_3[1].startswith('GGCTCCCCCC'))
        self.assertTrue(read_4_3[1].endswith('AGTTCTCGTG'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_4_1[2]), 1990)
            self.assertEqual(len(read_4_2[2]), 2868)
            self.assertEqual(len(read_4_3[2]), 1878)
            self.assertTrue(read_4_1[2].startswith("45.2*25-',"))
            self.assertTrue(read_4_1[2].endswith("$**+,,&$%3"))
            self.assertTrue(read_4_2[2].startswith("]213438::5"))
            self.assertTrue(read_4_2[2].endswith(".3/4:23]/]"))
            self.assertTrue(read_4_3[2].startswith("0430/6202,"))
            self.assertTrue(read_4_3[2].endswith("+2+,853322"))

        return read_type

    def test_results_normal_fastq(self):
        self.run_command('porechop -i INPUT -o OUTPUT.fastq')
        read_type = self.check_trimmed_reads()
        self.assertEqual(read_type, 'FASTQ')
        self.assertEqual(porechop.misc.get_compression_type(self.output_file), 'plain')

    def test_results_gz_fastq(self):
        self.run_command('porechop -i INPUT -o OUTPUT.fastq.gz')
        read_type = self.check_trimmed_reads()
        self.assertEqual(read_type, 'FASTQ')
        self.assertEqual(porechop.misc.get_compression_type(self.output_file), 'gz')

    def test_results_piped_fastq(self):
        self.run_command('porechop -i INPUT > OUTPUT.fastq')
        read_type = self.check_trimmed_reads()
        self.assertEqual(read_type, 'FASTQ')
        self.assertEqual(porechop.misc.get_compression_type(self.output_file), 'plain')

    def test_check_reads_1(self):
        """
        When only one read is checked, no adapters are found and nothing is trimmed.
        """
        out, _ = self.run_command('porechop -i INPUT -o OUTPUT.fastq --check_reads 1')
        self.assertTrue('No adapters found' in out)
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_1 = [x for x in trimmed_reads if x[0] == '1'][0]
        self.assertEqual(len(read_1[1]), 10000)
        read_2 = [x for x in trimmed_reads if x[0] == '2'][0]
        self.assertEqual(len(read_2[1]), 11000)
        read_3 = [x for x in trimmed_reads if x[0] == '3'][0]
        self.assertEqual(len(read_3[1]), 8000)

    def test_check_reads_2(self):
        """
        When only two reads are checked, only one adapter set is found.
        """
        out, _ = self.run_command('porechop -i INPUT -o OUTPUT.fastq --check_reads 2')
        self.assertTrue('Trimming adapters' in out)
        self.assertTrue('SQK-MAP006_Y_Top_SK63:' in out)
        self.assertTrue('SQK-MAP006_Y_Bottom_SK64:' in out)
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_1 = [x for x in trimmed_reads if x[0] == '1'][0]
        self.assertEqual(len(read_1[1]), 10000)
        read_2 = [x for x in trimmed_reads if x[0] == '2'][0]
        self.assertEqual(len(read_2[1]), 10970)
        read_3 = [x for x in trimmed_reads if x[0] == '3'][0]
        self.assertEqual(len(read_3[1]), 8000)

    def test_check_reads_3(self):
        """
        When three reads are checked, both adapter sets are found.
        """
        out, _ = self.run_command('porechop -i INPUT -o OUTPUT.fastq --check_reads 3')
        self.assertTrue('Trimming adapters' in out)
        self.assertTrue('SQK-MAP006_Y_Top_SK63:' in out)
        self.assertTrue('SQK-MAP006_Y_Bottom_SK64:' in out)
        self.assertTrue('SQK-NSK007_Y_Top:' in out)
        self.assertTrue('SQK-NSK007_Y_Bottom:' in out)
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_1 = [x for x in trimmed_reads if x[0] == '1'][0]
        self.assertEqual(len(read_1[1]), 10000)
        read_2 = [x for x in trimmed_reads if x[0] == '2'][0]
        self.assertEqual(len(read_2[1]), 10970)
        read_3 = [x for x in trimmed_reads if x[0] == '3'][0]
        self.assertEqual(len(read_3[1]), 7970)
