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


class TestOneAdapterSet(unittest.TestCase):
    """
    This test set contains only one adapter set: SQK-NSK007
    """
    def run_command(self, command):
        runner_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'porechop-runner.py')
        input_path = os.path.join(os.path.dirname(__file__), 'test_one_adapter_set.fastq')
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

    def test_verbosity_0_output(self):
        out, err = self.run_command('porechop -i INPUT -o OUTPUT.fastq --verbosity 0')
        self.assertEqual(out, '')
        self.assertEqual(err, '')

    def test_verbosity_1_output(self):
        out, err = self.run_command('porechop -i INPUT -o OUTPUT.fastq --verbosity 1')
        self.assertTrue('Trimming adapters from read ends' in out)
        self.assertTrue('SQK-NSK007_Y_Top:' in out)
        self.assertTrue('SQK-NSK007_Y_Bottom:' in out)
        self.assertTrue('4 / 9 reads' in out)
        self.assertTrue('3 / 9 reads' in out)
        self.assertEqual(err, '')

    def test_verbosity_2_output(self):
        out, err = self.run_command('porechop -i INPUT -o OUTPUT.fastq --verbosity 2')
        self.assertTrue('Trimming adapters from read ends' in out)
        self.assertTrue('4 / 9 reads' in out)
        self.assertTrue('3 / 9 reads' in out)
        self.assertTrue('CGCACCTCTCCCCTCTGCGTCCTAGGCACTAGATCCAAACCTAGTTCGCCTGAAATTTACTGATGCTAGACCG'
                        'AAACTTCGCGCCGACTACTCCATGGTT' in out)
        self.assertTrue('GCCCGTATCCACGTAAGAGTGCATCTCATTGCGCACAGGTATATCTGCCAGATAAGACGTCGAGG' in out)
        self.assertEqual(err, '')

    def test_piped_output(self):
        out, err = self.run_command('porechop -i INPUT')
        self.assertTrue('Trimming adapters from read ends' in err)
        self.assertTrue('SQK-NSK007_Y_Top:' in err)
        self.assertTrue('SQK-NSK007_Y_Bottom:' in err)
        self.assertTrue('4 / 9 reads' in err)
        self.assertTrue('3 / 9 reads' in err)
        self.assertTrue('@1\n' in out)
        self.assertTrue('@2\n' in out)
        self.assertTrue('@3\n' in out)

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
        self.assertTrue(read_1[1].startswith('CGCACCTCTC'))
        self.assertTrue(read_1[1].endswith('GATAGCAGGT'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_1[2]), 10000)
            self.assertTrue(read_1[2].startswith("4,)4+311.'"))
            self.assertTrue(read_1[2].endswith("3-16('...,"))

        # Read 2 has an adapter trimmed from the start.
        read_2 = trimmed_reads[1]
        self.assertEqual(read_2[0], '2')
        self.assertEqual(len(read_2[1]), 10970)
        self.assertTrue(read_2[1].startswith('GCATTGTGAG'))
        self.assertTrue(read_2[1].endswith('CAAGTGCCAG'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_2[2]), 10970)
            self.assertTrue(read_2[2].startswith("&)*-++$+,."))
            self.assertTrue(read_2[2].endswith("(*'%)%.2+)"))

        # Read 3 has an adapter trimmed from the end.
        read_3 = trimmed_reads[2]
        self.assertEqual(read_3[0], '3')
        self.assertEqual(len(read_3[1]), 7976)
        self.assertTrue(read_3[1].startswith('TTTTCGGCCG'))
        self.assertTrue(read_3[1].endswith('CGGCTGTTGG'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_3[2]), 7976)
            self.assertTrue(read_3[2].startswith("&-),1)5500"))
            self.assertTrue(read_3[2].endswith("14,)*,)+00"))

        # Read 4 has some adapter trimmed from both start and end.
        read_4 = trimmed_reads[3]
        self.assertEqual(read_4[0], '4')
        self.assertEqual(len(read_4[1]), 6968)
        self.assertTrue(read_4[1].startswith('TACCGTGTAA'))
        self.assertTrue(read_4[1].endswith('CAAGGGCAGC'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_4[2]), 6968)
            self.assertTrue(read_4[2].startswith("*3,+2&)'''"))
            self.assertTrue(read_4[2].endswith("+-)#'$$$'*"))

        # Read 5 has been split into two pieces based on an internal end adapter, but the second
        # half is discarded due to short length.
        read_5_1 = trimmed_reads[4]
        self.assertEqual(read_5_1[0], '5_1')
        self.assertEqual(len(read_5_1[1]), 3490)
        self.assertTrue(read_5_1[1].startswith('TCGGCTGTAC'))
        self.assertTrue(read_5_1[1].endswith('CGGATTGGGA'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_5_1[2]), 3490)
            self.assertTrue(read_5_1[2].startswith(")(1&172+(("))
            self.assertTrue(read_5_1[2].endswith(",/07389020"))

        # Read 6 has been split into two pieces based on an internal start adapter and both pieces
        # should remain.
        read_6_1 = trimmed_reads[5]
        read_6_2 = trimmed_reads[6]
        self.assertEqual(read_6_1[0], '6_1')
        self.assertEqual(read_6_2[0], '6_2')
        self.assertEqual(len(read_6_1[1]), 3900)
        self.assertEqual(len(read_6_2[1]), 7962)
        self.assertTrue(read_6_1[1].startswith('GTGGCCAATG'))
        self.assertTrue(read_6_1[1].endswith('ACTCAAGTTA'))
        self.assertTrue(read_6_2[1].startswith('AGGTTAGGTA'))
        self.assertTrue(read_6_2[1].endswith('GACGTGACGT'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_6_1[2]), 3900)
            self.assertEqual(len(read_6_2[2]), 7962)
            self.assertTrue(read_6_1[2].startswith("-.+,--.(/3"))
            self.assertTrue(read_6_1[2].endswith("2/(,-462,-"))
            self.assertTrue(read_6_2[2].startswith(")+)-/]++,("))
            self.assertTrue(read_6_2[2].endswith("33260,B(,-"))

        # Read 7 has a small piece of adapter trimmed from both start and end.
        read_7 = trimmed_reads[7]
        self.assertEqual(read_7[0], '7')
        self.assertEqual(len(read_7[1]), 1985)
        self.assertTrue(read_7[1].startswith('GAACTCGCGG'))
        self.assertTrue(read_7[1].endswith('TGGGACGCTA'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_7[2]), 1985)
            self.assertTrue(read_7[2].startswith("**(,6./7]."))
            self.assertTrue(read_7[2].endswith(")('*=]++15"))

        # Read 8 has been split into two pieces based on an internal start adapter as well as some
        # adapter at the start.
        read_8_1 = trimmed_reads[8]
        read_8_2 = trimmed_reads[9]
        self.assertEqual(read_8_1[0], '8_1')
        self.assertEqual(read_8_2[0], '8_2')
        self.assertEqual(len(read_8_1[1]), 4878)
        self.assertEqual(len(read_8_2[1]), 1962)
        self.assertTrue(read_8_1[1].startswith('ATAGTGAAGT'))
        self.assertTrue(read_8_1[1].endswith('GTCCTACTTA'))
        self.assertTrue(read_8_2[1].startswith('GACGTCAATA'))
        self.assertTrue(read_8_2[1].endswith('GTCTAAAACT'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_8_1[2]), 4878)
            self.assertEqual(len(read_8_2[2]), 1962)
            self.assertTrue(read_8_1[2].startswith(",,512;1151"))
            self.assertTrue(read_8_1[2].endswith("4%$+%%)'13"))
            self.assertTrue(read_8_2[2].startswith("].],'.*%+&"))
            self.assertTrue(read_8_2[2].endswith("]]020/1<:/"))

        # Read 9 has two internal adapters, but close together so the sequence in between is lost.
        read_9_1 = trimmed_reads[10]
        read_9_2 = trimmed_reads[11]
        self.assertEqual(read_9_1[0], '9_1')
        self.assertEqual(read_9_2[0], '9_2')
        self.assertEqual(len(read_9_1[1]), 1490)
        self.assertEqual(len(read_9_2[1]), 4262)
        self.assertTrue(read_9_1[1].startswith('CCTTGACAAG'))
        self.assertTrue(read_9_1[1].endswith('CAGCGCGTAT'))
        self.assertTrue(read_9_2[1].startswith('CAGGTTGTGT'))
        self.assertTrue(read_9_2[1].endswith('GGATACTACC'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_9_1[2]), 1490)
            self.assertEqual(len(read_9_2[2]), 4262)
            self.assertTrue(read_9_1[2].startswith(".0F]31557-"))
            self.assertTrue(read_9_1[2].endswith("4&(+343-'0"))
            self.assertTrue(read_9_2[2].startswith("$%*218.2*0"))
            self.assertTrue(read_9_2[2].endswith(".,//()&(&*"))

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

    def test_results_normal_fasta(self):
        self.run_command('porechop -i INPUT -o OUTPUT.fasta')
        read_type = self.check_trimmed_reads()
        self.assertEqual(read_type, 'FASTA')
        self.assertEqual(porechop.misc.get_compression_type(self.output_file), 'plain')

    def test_results_gz_fasta(self):
        self.run_command('porechop -i INPUT -o OUTPUT.fasta.gz')
        read_type = self.check_trimmed_reads()
        self.assertEqual(read_type, 'FASTA')
        self.assertEqual(porechop.misc.get_compression_type(self.output_file), 'gz')

    def test_results_piped_fasta(self):
        self.run_command('porechop -i INPUT --format fasta > OUTPUT.fasta')
        read_type = self.check_trimmed_reads()
        self.assertEqual(read_type, 'FASTA')
        self.assertEqual(porechop.misc.get_compression_type(self.output_file), 'plain')

    def test_single_threaded(self):
        self.run_command('porechop -i INPUT -o OUTPUT.fastq --threads 1')
        self.check_trimmed_reads()

    def test_multi_threaded(self):
        self.run_command('porechop -i INPUT -o OUTPUT.fastq --threads 8')
        self.check_trimmed_reads()

    def test_end_size_1(self):
        self.run_command('porechop -i INPUT -o OUTPUT.fastq --end_size 50')
        self.check_trimmed_reads()

    def test_end_size_2(self):
        self.run_command('porechop -i INPUT -o OUTPUT.fastq --end_size 100')
        self.check_trimmed_reads()

    def test_end_size_3(self):
        self.run_command('porechop -i INPUT -o OUTPUT.fastq --end_size 200')
        self.check_trimmed_reads()

    def test_min_trim_size_1(self):
        """
        Increasing this setting will filter out the adapter trimming from read 7.
        """
        self.run_command('porechop -i INPUT -o OUTPUT.fastq --min_trim_size 5')
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_7 = [x for x in trimmed_reads if x[0] == '7'][0]
        self.assertEqual(read_7[0], '7')
        self.assertEqual(len(read_7[1]), 1985)
        self.assertTrue(read_7[1].startswith('GAACTCGCGG'))
        self.assertTrue(read_7[1].endswith('TGGGACGCTA'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_7[2]), 1985)
            self.assertTrue(read_7[2].startswith("**(,6./7]."))
            self.assertTrue(read_7[2].endswith(")('*=]++15"))

    def test_min_trim_size_2(self):
        self.run_command('porechop -i INPUT -o OUTPUT.fastq --min_trim_size 6')
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_7 = [x for x in trimmed_reads if x[0] == '7'][0]
        self.assertEqual(read_7[0], '7')
        self.assertEqual(len(read_7[1]), 1992)
        self.assertTrue(read_7[1].startswith('GAACTCGCGG'))
        self.assertTrue(read_7[1].endswith('CTACCGCAAT'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_7[2]), 1992)
            self.assertTrue(read_7[2].startswith("**(,6./7]."))
            self.assertTrue(read_7[2].endswith("+1505'+,),"))

    def test_min_trim_size_3(self):
        self.run_command('porechop -i INPUT -o OUTPUT.fastq --min_trim_size 7')
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_7 = [x for x in trimmed_reads if x[0] == '7'][0]
        self.assertEqual(read_7[0], '7')
        self.assertEqual(len(read_7[1]), 2000)
        self.assertTrue(read_7[1].startswith('ATTGCTAGGA'))
        self.assertTrue(read_7[1].endswith('CTACCGCAAT'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_7[2]), 2000)
            self.assertTrue(read_7[2].startswith("/-]%+(-***"))
            self.assertTrue(read_7[2].endswith("+1505'+,),"))

    def test_extra_middle_trim_1(self):
        self.run_command('porechop -i INPUT -o OUTPUT.fastq --extra_middle_trim_good_side 0'
                         ' --extra_middle_trim_bad_side 0')
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_6_1 = [x for x in trimmed_reads if x[0] == '6_1'][0]
        read_6_2 = [x for x in trimmed_reads if x[0] == '6_2'][0]
        self.assertEqual(len(read_6_1[1]), 4000)
        self.assertEqual(len(read_6_2[1]), 7972)

    def test_extra_middle_trim_2(self):
        self.run_command('porechop -i INPUT -o OUTPUT.fastq --extra_middle_trim_good_side 123'
                         ' --extra_middle_trim_bad_side 456')
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_6_1 = [x for x in trimmed_reads if x[0] == '6_1'][0]
        read_6_2 = [x for x in trimmed_reads if x[0] == '6_2'][0]
        self.assertEqual(len(read_6_1[1]), 3544)
        self.assertEqual(len(read_6_2[1]), 7849)

    def test_middle_threshold_1(self):
        """
        Read 6 has an inexact match for the start adapter, so increasing the middle_threshold
        option should prevent the split.
        """
        self.run_command('porechop -i INPUT -o OUTPUT.fastq --middle_threshold 96')
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_6_1 = [x for x in trimmed_reads if x[0] == '6_1'][0]
        read_6_2 = [x for x in trimmed_reads if x[0] == '6_2'][0]
        self.assertEqual(len(read_6_1[1]), 3900)
        self.assertEqual(len(read_6_2[1]), 7962)

    def test_middle_threshold_2(self):
        self.run_command('porechop -i INPUT -o OUTPUT.fastq --middle_threshold 97')
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_6 = [x for x in trimmed_reads if x[0] == '6'][0]
        self.assertEqual(len(read_6[1]), 12000)

    def test_check_reads(self):
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
        read_4 = [x for x in trimmed_reads if x[0] == '4'][0]
        self.assertEqual(len(read_4[1]), 7000)
        read_5 = [x for x in trimmed_reads if x[0] == '5'][0]
        self.assertEqual(len(read_5[1]), 4000)
        read_6 = [x for x in trimmed_reads if x[0] == '6'][0]
        self.assertEqual(len(read_6[1]), 12000)
        read_7 = [x for x in trimmed_reads if x[0] == '7'][0]
        self.assertEqual(len(read_7[1]), 2000)
        read_8 = [x for x in trimmed_reads if x[0] == '8'][0]
        self.assertEqual(len(read_8[1]), 7000)
        read_9 = [x for x in trimmed_reads if x[0] == '9'][0]
        self.assertEqual(len(read_9[1]), 6000)

    def test_adapter_threshold_1(self):
        """
        When only two reads are checked with an adapter threshold of 100%, no adapters are found
        (because the adapter in read 2 is not an exact match).
        """
        out, _ = self.run_command('porechop -i INPUT -o OUTPUT.fastq --check_reads 2 '
                                  '--adapter_threshold 100')
        self.assertTrue('No adapters found' in out)
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_2 = [x for x in trimmed_reads if x[0] == '2'][0]
        self.assertEqual(len(read_2[1]), 11000)

    def test_adapter_threshold_2(self):
        """
        When the adapter threshold is lowered to 90%, Porechop finds the adapter in read 2..
        """
        out, _ = self.run_command('porechop -i INPUT -o OUTPUT.fastq --check_reads 2 '
                                  '--adapter_threshold 90')
        self.assertTrue('No adapters found' not in out)
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_2 = [x for x in trimmed_reads if x[0] == '2'][0]
        self.assertEqual(len(read_2[1]), 10970)

    def test_extra_end_trim_1(self):
        out, _ = self.run_command('porechop -i INPUT -o OUTPUT.fastq --extra_end_trim 0')
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_2 = [x for x in trimmed_reads if x[0] == '2'][0]
        self.assertEqual(len(read_2[1]), 10972)
        read_3 = [x for x in trimmed_reads if x[0] == '3'][0]
        self.assertEqual(len(read_3[1]), 7978)
        read_7 = [x for x in trimmed_reads if x[0] == '7'][0]
        self.assertEqual(len(read_7[1]), 1989)

    def test_extra_end_trim_2(self):
        out, _ = self.run_command('porechop -i INPUT -o OUTPUT.fastq --extra_end_trim 10')
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_2 = [x for x in trimmed_reads if x[0] == '2'][0]
        self.assertEqual(len(read_2[1]), 10962)
        read_3 = [x for x in trimmed_reads if x[0] == '3'][0]
        self.assertEqual(len(read_3[1]), 7968)
        read_7 = [x for x in trimmed_reads if x[0] == '7'][0]
        self.assertEqual(len(read_7[1]), 1969)

    def test_extra_end_trim_3(self):
        out, _ = self.run_command('porechop -i INPUT -o OUTPUT.fastq --extra_end_trim 100')
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_2 = [x for x in trimmed_reads if x[0] == '2'][0]
        self.assertEqual(len(read_2[1]), 10872)
        read_3 = [x for x in trimmed_reads if x[0] == '3'][0]
        self.assertEqual(len(read_3[1]), 7878)
        read_7 = [x for x in trimmed_reads if x[0] == '7'][0]
        self.assertEqual(len(read_7[1]), 1789)

    def test_min_split_read_size_1(self):
        out, _ = self.run_command('porechop -i INPUT -o OUTPUT.fastq --min_split_read_size 1 '
                                  '--extra_middle_trim_good_side 0 --extra_middle_trim_bad_side 0')
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_5_1 = [x for x in trimmed_reads if x[0] == '5_1'][0]
        read_5_2 = [x for x in trimmed_reads if x[0] == '5_2'][0]
        self.assertEqual(len(read_5_1[1]), 3500)
        self.assertEqual(len(read_5_2[1]), 478)

    def test_min_split_read_size_2(self):
        out, _ = self.run_command('porechop -i INPUT -o OUTPUT.fastq --min_split_read_size 478 '
                                  '--extra_middle_trim_good_side 0 --extra_middle_trim_bad_side 0')
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_5_1 = [x for x in trimmed_reads if x[0] == '5_1'][0]
        read_5_2 = [x for x in trimmed_reads if x[0] == '5_2'][0]
        self.assertEqual(len(read_5_1[1]), 3500)
        self.assertEqual(len(read_5_2[1]), 478)

    def test_min_split_read_size_3(self):
        out, _ = self.run_command('porechop -i INPUT -o OUTPUT.fastq --min_split_read_size 479 '
                                  '--extra_middle_trim_good_side 0 --extra_middle_trim_bad_side 0')
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_5_1 = [x for x in trimmed_reads if x[0] == '5_1'][0]
        self.assertEqual(len(read_5_1[1]), 3500)
        self.assertEqual(len([x for x in trimmed_reads if x[0] == '5_2']), 0)

    def test_min_split_read_size_4(self):
        out, _ = self.run_command('porechop -i INPUT -o OUTPUT.fastq --min_split_read_size 1 '
                                  '--extra_middle_trim_good_side 50 '
                                  '--extra_middle_trim_bad_side 50')
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_5_1 = [x for x in trimmed_reads if x[0] == '5_1'][0]
        read_5_2 = [x for x in trimmed_reads if x[0] == '5_2'][0]
        self.assertEqual(len(read_5_1[1]), 3450)
        self.assertEqual(len(read_5_2[1]), 428)

    def test_min_split_read_size_5(self):
        out, _ = self.run_command('porechop -i INPUT -o OUTPUT.fastq --min_split_read_size 428 '
                                  '--extra_middle_trim_good_side 50 '
                                  '--extra_middle_trim_bad_side 50')
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_5_1 = [x for x in trimmed_reads if x[0] == '5_1'][0]
        read_5_2 = [x for x in trimmed_reads if x[0] == '5_2'][0]
        self.assertEqual(len(read_5_1[1]), 3450)
        self.assertEqual(len(read_5_2[1]), 428)

    def test_min_split_read_size_6(self):
        out, _ = self.run_command('porechop -i INPUT -o OUTPUT.fastq --min_split_read_size 429 '
                                  '--extra_middle_trim_good_side 50 '
                                  '--extra_middle_trim_bad_side 50')
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_5_1 = [x for x in trimmed_reads if x[0] == '5_1'][0]
        self.assertEqual(len(read_5_1[1]), 3450)
        self.assertEqual(len([x for x in trimmed_reads if x[0] == '5_2']), 0)

    def test_min_split_read_size_7(self):
        out, _ = self.run_command('porechop -i INPUT -o OUTPUT.fastq --min_split_read_size 1 '
                                  '--extra_middle_trim_good_side 0 --extra_middle_trim_bad_side 0')
        trimmed_reads, read_type = self.load_trimmed_reads()
        read_9_1 = [x for x in trimmed_reads if x[0] == '9_1'][0]
        read_9_2 = [x for x in trimmed_reads if x[0] == '9_2'][0]
        read_9_3 = [x for x in trimmed_reads if x[0] == '9_3'][0]
        self.assertEqual(len(read_9_1[1]), 1500)
        self.assertEqual(len(read_9_2[1]), 178)
        self.assertEqual(len(read_9_3[1]), 4272)
