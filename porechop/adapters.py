"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Porechop

This module contains the class and sequences for known adapters used in Oxford Nanopore library
preparation kits.

This file is part of Porechop. Porechop is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Porechop is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Porechop. If
not, see <http://www.gnu.org/licenses/>.
"""


class Adapter(object):

    def __init__(self, name, start_sequence=None, end_sequence=None, both_ends_sequence=None):
        self.name = name
        self.start_sequence = start_sequence if start_sequence else []
        self.end_sequence = end_sequence if end_sequence else []
        if both_ends_sequence:
            self.start_sequence = both_ends_sequence
            self.end_sequence = both_ends_sequence
        self.best_start_score, self.best_end_score = 0.0, 0.0

    def best_start_or_end_score(self):
        return max(self.best_start_score, self.best_end_score)


ADAPTERS = [Adapter('SQK-MAP006',
                    start_sequence=('SQK-MAP006_Y_Top_SK63',    'GGTTGTTTCTGTTGGTGCTGATATTGCT'),
                    end_sequence=  ('SQK-MAP006_Y_Bottom_SK64', 'GCAATATCAGCACCAACAGAAA')),
            Adapter('SQK-MAP006 Short',
                    start_sequence=('SQK-MAP006_Short_Y_Top_LI32',    'CGGCGTCTGCTTGGGTGTTTAACCT'),
                    end_sequence=  ('SQK-MAP006_Short_Y_Bottom_LI33', 'GGTTAAACACCCAAGCAGACGCCG')),
            Adapter('SQK-NSK007',
                    start_sequence=('SQK-NSK007_Y_Top',    'AATGTACTTCGTTCAGTTACGTATTGCT'),
                    end_sequence=  ('SQK-NSK007_Y_Bottom', 'GCAATACGTAACTGAACGAAGT')),
            Adapter('Native barcoding 1',
                    start_sequence=('NB01_rev', 'AGGTTAACACAAAGACACCGACAACTTTCTTCAGCACC'),
                    end_sequence=  ('NB01',     'GGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCT')),
            Adapter('Native barcoding 2',
                    start_sequence=('NB02_rev', 'AGGTTAAACAGACGACTACAAACGGAATCGACAGCACC'),
                    end_sequence=  ('NB02',     'GGTGCTGTCGATTCCGTTTGTAGTCGTCTGTTTAACCT')),
            Adapter('Native barcoding 3',
                    start_sequence=('NB03_rev', 'AGGTTAACCTGGTAACTGGGACACAAGACTCCAGCACC'),
                    end_sequence=  ('NB03',     'GGTGCTGGAGTCTTGTGTCCCAGTTACCAGGTTAACCT')),
            Adapter('Native barcoding 4',
                    start_sequence=('NB04_rev', 'AGGTTAATAGGGAAACACGATAGAATCCGAACAGCACC'),
                    end_sequence=  ('NB04',     'GGTGCTGTTCGGATTCTATCGTGTTTCCCTATTAACCT')),
            Adapter('Native barcoding 5',
                    start_sequence=('NB05_rev', 'AGGTTAAAAGGTTACACAAACCCTGGACAAGCAGCACC'),
                    end_sequence=  ('NB05',     'GGTGCTGCTTGTCCAGGGTTTGTGTAACCTTTTAACCT')),
            Adapter('Native barcoding 6',
                    start_sequence=('NB06_rev', 'AGGTTAAGACTACTTTCTGCCTTTGCGAGAACAGCACC'),
                    end_sequence=  ('NB06',     'GGTGCTGTTCTCGCAAAGGCAGAAAGTAGTCTTAACCT')),
            Adapter('Native barcoding 7',
                    start_sequence=('NB07_rev', 'AGGTTAAAAGGATTCATTCCCACGGTAACACCAGCACC'),
                    end_sequence=  ('NB07',     'GGTGCTGGTGTTACCGTGGGAATGAATCCTTTTAACCT')),
            Adapter('Native barcoding 8',
                    start_sequence=('NB08_rev', 'AGGTTAAACGTAACTTGGTTTGTTCCCTGAACAGCACC'),
                    end_sequence=  ('NB08',     'GGTGCTGTTCAGGGAACAAACCAAGTTACGTTTAACCT')),
            Adapter('Native barcoding 9',
                    start_sequence=('NB09_rev', 'AGGTTAAAACCAAGACTCGCTGTGCCTAGTTCAGCACC'),
                    end_sequence=  ('NB09',     'GGTGCTGAACTAGGCACAGCGAGTCTTGGTTTTAACCT')),
            Adapter('Native barcoding 10',
                    start_sequence=('NB10_rev', 'AGGTTAAGAGAGGACAAAGGTTTCAACGCTTCAGCACC'),
                    end_sequence=  ('NB10',     'GGTGCTGAAGCGTTGAAACCTTTGTCCTCTCTTAACCT')),
            Adapter('Native barcoding 11',
                    start_sequence=('NB11_rev', 'AGGTTAATCCATTCCCTCCGATAGATGAAACCAGCACC'),
                    end_sequence=  ('NB11',     'GGTGCTGGTTTCATCTATCGGAGGGAATGGATTAACCT')),
            Adapter('Native barcoding 12',
                    start_sequence=('NB12_rev', 'AGGTTAATCCGATTCTGCTTCTTTCTACCTGCAGCACC'),
                    end_sequence=  ('NB12',     'GGTGCTGCAGGTAGAAAGAAGCAGAATCGGATTAACCT')),
            ]
