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


class TestAlbacoreDirectory(unittest.TestCase):
    """
    Tests barcode demultiplexing from an Albacore output directory.
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

    def test_albacore_directory_1(self):
        """
        When just the barcode01 reads are given as a single file, Porechop doesn't treat it as
        Albacore output.
        """
        out, _ = self.run_command('porechop -i INPUT -b BARCODE_DIR',
                                  'test_albacore_directory/workspace/barcode01/' +
                                  'fastq_runid_e1df3d311cc267665db81df4aa863d073b0494ae_0.fastq')

        bc01_trimmed_reads = self.load_trimmed_reads('BC01.fastq')
        bc03_trimmed_reads = self.load_trimmed_reads('BC03.fastq')
        none_trimmed_reads = self.load_trimmed_reads('none.fastq')

        bc01_read_names = sorted(x[0] for x in bc01_trimmed_reads)
        bc03_read_names = sorted(x[0] for x in bc03_trimmed_reads)
        none_read_names = sorted(x[0] for x in none_trimmed_reads)

        self.assertEqual(bc01_read_names, ['432342b4-d6fb-4915-96ca-1b4a570779cf',
                                           'a4a5afa1-747e-47ea-abb9-c54029c315ec',
                                           'd82a6ff9-99cf-4f23-9246-cedfc6e8c88b'])
        self.assertEqual(bc03_read_names, ['b8db6753-a98a-40c3-8f80-aadedcf4d1e4'])
        self.assertEqual(none_read_names, ['04a00c9c-4a3a-4ba6-8161-ba01d9a77f22',
                                           '9dc8c47d-e574-497b-b24d-5cf8f09cb189',
                                           'ba8d33e4-7b53-42e8-8e6f-1d017698a21e',
                                           'f38d7656-735d-4be9-92ce-94a568bada56',
                                           'ffa40714-6027-4ebc-9f49-ecde4f6511f3'])

    def test_albacore_directory_2(self):
        """
        When just the barcode02 reads are given as a single file, Porechop doesn't treat it as
        Albacore output.
        """
        out, _ = self.run_command('porechop -i INPUT -b BARCODE_DIR',
                                  'test_albacore_directory/workspace/barcode02/' +
                                  'fastq_runid_e1df3d311cc267665db81df4aa863d073b0494ae_0.fastq')

        bc02_trimmed_reads = self.load_trimmed_reads('BC02.fastq')
        bc03_trimmed_reads = self.load_trimmed_reads('BC03.fastq')
        none_trimmed_reads = self.load_trimmed_reads('none.fastq')

        bc02_read_names = sorted(x[0] for x in bc02_trimmed_reads)
        bc03_read_names = sorted(x[0] for x in bc03_trimmed_reads)
        none_read_names = sorted(x[0] for x in none_trimmed_reads)

        self.assertEqual(bc02_read_names, ['15c8df94-6795-451b-8f71-bc32580068f2',
                                           '28e41a6a-25db-4023-a27e-f3ffb59f238a',
                                           '4dd6464b-5cea-4f3a-8f26-740983ca2a57'])
        self.assertEqual(bc03_read_names, ['ca5efc50-1e67-4df4-95e9-afe4a5dbf35c'])
        self.assertEqual(none_read_names, ['1c25347a-7a36-4f3a-a5a2-5339ae585d4a',
                                           '5b0c7ed8-c292-4948-a235-c3644a15855b',
                                           'de7b3819-a55c-4c6c-88ea-383de089ebab'])

    def test_albacore_directory_3(self):
        """
        When just the barcode03 reads are given as a single file, Porechop doesn't treat it as
        Albacore output.
        """
        out, _ = self.run_command('porechop -i INPUT -b BARCODE_DIR',
                                  'test_albacore_directory/workspace/barcode03/' +
                                  'fastq_runid_e1df3d311cc267665db81df4aa863d073b0494ae_0.fastq')

        bc01_trimmed_reads = self.load_trimmed_reads('BC01.fastq')
        bc03_trimmed_reads = self.load_trimmed_reads('BC03.fastq')
        none_trimmed_reads = self.load_trimmed_reads('none.fastq')

        bc01_read_names = sorted(x[0] for x in bc01_trimmed_reads)
        bc03_read_names = sorted(x[0] for x in bc03_trimmed_reads)
        none_read_names = sorted(x[0] for x in none_trimmed_reads)

        self.assertEqual(bc01_read_names, ['283ab5e9-492e-4b28-a200-0362b4755ced'])
        self.assertEqual(bc03_read_names, ['47277a55-611f-415e-b0d4-3d3326e0e125',
                                           '73d0e477-fbfe-49c2-abc9-4974dd765bc8',
                                           'd75205b6-be37-4a5c-bed3-d6b218e7579b'])
        self.assertEqual(none_read_names, ['618187f9-2f66-4bf0-a1b6-e735b7024e86',
                                           '6f612e35-8174-4dc9-a239-503b10915267',
                                           'd30719ed-c801-48db-939c-98b10137b414',
                                           'fb5162ca-4645-4a0d-9336-3aa5fdf18d74'])

    def test_albacore_directory_unclassified(self):
        """
        When just the unclassified reads are given as a single file, Porechop doesn't treat it as
        Albacore output.
        """
        out, _ = self.run_command('porechop -i INPUT -b BARCODE_DIR',
                                  'test_albacore_directory/workspace/unclassified/' +
                                  'fastq_runid_e1df3d311cc267665db81df4aa863d073b0494ae_0.fastq')

        bc01_trimmed_reads = self.load_trimmed_reads('BC01.fastq')
        bc02_trimmed_reads = self.load_trimmed_reads('BC02.fastq')
        bc03_trimmed_reads = self.load_trimmed_reads('BC03.fastq')
        none_trimmed_reads = self.load_trimmed_reads('none.fastq')

        bc01_read_names = sorted(x[0] for x in bc01_trimmed_reads)
        bc02_read_names = sorted(x[0] for x in bc02_trimmed_reads)
        bc03_read_names = sorted(x[0] for x in bc03_trimmed_reads)
        none_read_names = sorted(x[0] for x in none_trimmed_reads)

        self.assertEqual(bc01_read_names, ['a7ffe707-2733-43b9-a9df-950566189723'])
        self.assertEqual(bc02_read_names, ['111fc53c-1941-4c5b-82e2-c80dc0962bd7',
                                           '50e78d16-e5f9-4b77-8ce4-dde6fc6eea39'])
        self.assertEqual(bc03_read_names, ['7557893f-a3a1-4294-96b7-89f9e9d4470c'])
        self.assertEqual(none_read_names, ['4d1ff7d1-b2f1-4ecc-8990-09edf0aa1970',
                                           '5e656a69-bd50-45c9-8177-a2aa38d0b04f',
                                           '7680810e-db1f-4cde-b442-642d7ac024ca'])

    def test_albacore_directory_all(self):
        """
        When the whole directory is given as input, it's treated as an Albacore directory and reads
        where Porechop and Albacore disagree are put into 'none'.
        """
        out, _ = self.run_command('porechop -i INPUT -b BARCODE_DIR', 'test_albacore_directory')

        bc01_trimmed_reads = self.load_trimmed_reads('BC01.fastq')
        bc02_trimmed_reads = self.load_trimmed_reads('BC02.fastq')
        bc03_trimmed_reads = self.load_trimmed_reads('BC03.fastq')
        none_trimmed_reads = self.load_trimmed_reads('none.fastq')

        bc01_read_names = sorted(x[0] for x in bc01_trimmed_reads)
        bc02_read_names = sorted(x[0] for x in bc02_trimmed_reads)
        bc03_read_names = sorted(x[0] for x in bc03_trimmed_reads)
        none_read_names = sorted(x[0] for x in none_trimmed_reads)

        self.assertEqual(bc01_read_names, ['432342b4-d6fb-4915-96ca-1b4a570779cf',
                                           'a4a5afa1-747e-47ea-abb9-c54029c315ec',
                                           'd82a6ff9-99cf-4f23-9246-cedfc6e8c88b'])
        self.assertEqual(bc02_read_names, ['15c8df94-6795-451b-8f71-bc32580068f2',
                                           '28e41a6a-25db-4023-a27e-f3ffb59f238a',
                                           '4dd6464b-5cea-4f3a-8f26-740983ca2a57'])
        self.assertEqual(bc03_read_names, ['47277a55-611f-415e-b0d4-3d3326e0e125',
                                           '73d0e477-fbfe-49c2-abc9-4974dd765bc8',
                                           'd75205b6-be37-4a5c-bed3-d6b218e7579b'])
        self.assertEqual(none_read_names, ['04a00c9c-4a3a-4ba6-8161-ba01d9a77f22',
                                           '111fc53c-1941-4c5b-82e2-c80dc0962bd7',
                                           '1c25347a-7a36-4f3a-a5a2-5339ae585d4a',
                                           '283ab5e9-492e-4b28-a200-0362b4755ced',
                                           '4d1ff7d1-b2f1-4ecc-8990-09edf0aa1970',
                                           '50e78d16-e5f9-4b77-8ce4-dde6fc6eea39',
                                           '5b0c7ed8-c292-4948-a235-c3644a15855b',
                                           '5e656a69-bd50-45c9-8177-a2aa38d0b04f',
                                           '618187f9-2f66-4bf0-a1b6-e735b7024e86',
                                           '6f612e35-8174-4dc9-a239-503b10915267',
                                           '7557893f-a3a1-4294-96b7-89f9e9d4470c',
                                           '7680810e-db1f-4cde-b442-642d7ac024ca',
                                           '9dc8c47d-e574-497b-b24d-5cf8f09cb189',
                                           'a7ffe707-2733-43b9-a9df-950566189723',
                                           'b8db6753-a98a-40c3-8f80-aadedcf4d1e4',
                                           'ba8d33e4-7b53-42e8-8e6f-1d017698a21e',
                                           'ca5efc50-1e67-4df4-95e9-afe4a5dbf35c',
                                           'd30719ed-c801-48db-939c-98b10137b414',
                                           'de7b3819-a55c-4c6c-88ea-383de089ebab',
                                           'f38d7656-735d-4be9-92ce-94a568bada56',
                                           'fb5162ca-4645-4a0d-9336-3aa5fdf18d74',
                                           'ffa40714-6027-4ebc-9f49-ecde4f6511f3'])
