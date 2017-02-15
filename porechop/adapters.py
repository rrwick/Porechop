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
                    start_sequence=('NB01',     'GGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCT'),
                    end_sequence=  ('NB01_rev', 'AGGTTAACACAAAGACACCGACAACTTTCTTCAGCACC')),
            Adapter('Native barcoding 2',
                    start_sequence=('NB02',     'GGTGCTGTCGATTCCGTTTGTAGTCGTCTGTTTAACCT'),
                    end_sequence=  ('NB02_rev', 'AGGTTAAACAGACGACTACAAACGGAATCGACAGCACC')),
            Adapter('Native barcoding 3',
                    start_sequence=('NB03',     'GGTGCTGGAGTCTTGTGTCCCAGTTACCAGGTTAACCT'),
                    end_sequence=  ('NB03_rev', 'AGGTTAACCTGGTAACTGGGACACAAGACTCCAGCACC')),
            Adapter('Native barcoding 4',
                    start_sequence=('NB04',     'GGTGCTGTTCGGATTCTATCGTGTTTCCCTATTAACCT'),
                    end_sequence=  ('NB04_rev', 'AGGTTAATAGGGAAACACGATAGAATCCGAACAGCACC')),
            Adapter('Native barcoding 5',
                    start_sequence=('NB05',     'GGTGCTGCTTGTCCAGGGTTTGTGTAACCTTTTAACCT'),
                    end_sequence=  ('NB05_rev', 'AGGTTAAAAGGTTACACAAACCCTGGACAAGCAGCACC')),
            Adapter('Native barcoding 6',
                    start_sequence=('NB06',     'GGTGCTGTTCTCGCAAAGGCAGAAAGTAGTCTTAACCT'),
                    end_sequence=  ('NB06_rev', 'AGGTTAAGACTACTTTCTGCCTTTGCGAGAACAGCACC')),
            Adapter('Native barcoding 7',
                    start_sequence=('NB07',     'GGTGCTGGTGTTACCGTGGGAATGAATCCTTTTAACCT'),
                    end_sequence=  ('NB07_rev', 'AGGTTAAAAGGATTCATTCCCACGGTAACACCAGCACC')),
            Adapter('Native barcoding 8',
                    start_sequence=('NB08',     'GGTGCTGTTCAGGGAACAAACCAAGTTACGTTTAACCT'),
                    end_sequence=  ('NB08_rev', 'AGGTTAAACGTAACTTGGTTTGTTCCCTGAACAGCACC')),
            Adapter('Native barcoding 9',
                    start_sequence=('NB09',     'GGTGCTGAACTAGGCACAGCGAGTCTTGGTTTTAACCT'),
                    end_sequence=  ('NB09_rev', 'AGGTTAAAACCAAGACTCGCTGTGCCTAGTTCAGCACC')),
            Adapter('Native barcoding 10',
                    start_sequence=('NB10',     'GGTGCTGAAGCGTTGAAACCTTTGTCCTCTCTTAACCT'),
                    end_sequence=  ('NB10_rev', 'AGGTTAAGAGAGGACAAAGGTTTCAACGCTTCAGCACC')),
            Adapter('Native barcoding 11',
                    start_sequence=('NB11',     'GGTGCTGGTTTCATCTATCGGAGGGAATGGATTAACCT'),
                    end_sequence=  ('NB11_rev', 'AGGTTAATCCATTCCCTCCGATAGATGAAACCAGCACC')),
            Adapter('Native barcoding 12',
                    start_sequence=('NB12',     'GGTGCTGCAGGTAGAAAGAAGCAGAATCGGATTAACCT'),
                    end_sequence=  ('NB12_rev', 'AGGTTAATCCGATTCTGCTTCTTTCTACCTGCAGCACC')),
            Adapter('PCR barcoding 1',
                    start_sequence=('BC01',     'AAGAAAGTTGTCGGTGTCTTTGTG'),
                    end_sequence=  ('BC01_rev', 'CACAAAGACACCGACAACTTTCTT')),
            Adapter('PCR barcoding 2',
                    start_sequence=('BC02',     'TCGATTCCGTTTGTAGTCGTCTGT'),
                    end_sequence=  ('BC02_rev', 'ACAGACGACTACAAACGGAATCGA')),
            Adapter('PCR barcoding 3',
                    start_sequence=('BC03',     'GAGTCTTGTGTCCCAGTTACCAGG'),
                    end_sequence=  ('BC03_rev', 'CCTGGTAACTGGGACACAAGACTC')),
            Adapter('PCR barcoding 4',
                    start_sequence=('BC04',     'TTCGGATTCTATCGTGTTTCCCTA'),
                    end_sequence=  ('BC04_rev', 'TAGGGAAACACGATAGAATCCGAA')),
            Adapter('PCR barcoding 5',
                    start_sequence=('BC05',     'CTTGTCCAGGGTTTGTGTAACCTT'),
                    end_sequence=  ('BC05_rev', 'AAGGTTACACAAACCCTGGACAAG')),
            Adapter('PCR barcoding 6',
                    start_sequence=('BC06',     'TTCTCGCAAAGGCAGAAAGTAGTC'),
                    end_sequence=  ('BC06_rev', 'GACTACTTTCTGCCTTTGCGAGAA')),
            Adapter('PCR barcoding 7',
                    start_sequence=('BC07',     'GTGTTACCGTGGGAATGAATCCTT'),
                    end_sequence=  ('BC07_rev', 'AAGGATTCATTCCCACGGTAACAC')),
            Adapter('PCR barcoding 8',
                    start_sequence=('BC08',     'TTCAGGGAACAAACCAAGTTACGT'),
                    end_sequence=  ('BC08_rev', 'ACGTAACTTGGTTTGTTCCCTGAA')),
            Adapter('PCR barcoding 9',
                    start_sequence=('BC09',     'AACTAGGCACAGCGAGTCTTGGTT'),
                    end_sequence=  ('BC09_rev', 'AACCAAGACTCGCTGTGCCTAGTT')),
            Adapter('PCR barcoding 10',
                    start_sequence=('BC10',     'AAGCGTTGAAACCTTTGTCCTCTC'),
                    end_sequence=  ('BC10_rev', 'GAGAGGACAAAGGTTTCAACGCTT')),
            Adapter('PCR barcoding 11',
                    start_sequence=('BC11',     'GTTTCATCTATCGGAGGGAATGGA'),
                    end_sequence=  ('BC11_rev', 'TCCATTCCCTCCGATAGATGAAAC')),
            Adapter('PCR barcoding 12',
                    start_sequence=('BC12',     'CAGGTAGAAAGAAGCAGAATCGGA'),
                    end_sequence=  ('BC12_rev', 'TCCGATTCTGCTTCTTTCTACCTG'))
            ]
