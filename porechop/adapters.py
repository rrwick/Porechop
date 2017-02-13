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
                    start_sequence=('SQK-MAP006_Y_Top_SK63', 'GGTTGTTTCTGTTGGTGCTGATATTGCT'),
                    end_sequence=('SQK-MAP006_Y_Bottom_SK64', 'GCAATATCAGCACCAACAGAAA')),
            Adapter('SQK-MAP006 Short',
                    start_sequence=('SQK-MAP006_Short_Y_Top_LI32', 'CGGCGTCTGCTTGGGTGTTTAACCT'),
                    end_sequence=('SQK-MAP006_Short_Y_Bottom_LI33', 'GGTTAAACACCCAAGCAGACGCCG')),
            Adapter('SQK-NSK007',
                    start_sequence=('SQK-NSK007_Y_Top', 'AATGTACTTCGTTCAGTTACGTATTGCT'),
                    end_sequence=('SQK-NSK007_Y_Bottom', 'GCAATACGTAACTGAACGAAGT')),
            Adapter('Native barcoding 1',
                    both_ends_sequence=('NB01', 'GGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCT')),
            Adapter('Native barcoding 2',
                    both_ends_sequence=('NB02', 'GGTGCTGTCGATTCCGTTTGTAGTCGTCTGTTTAACCT')),
            Adapter('Native barcoding 3',
                    both_ends_sequence=('NB03', 'GGTGCTGGAGTCTTGTGTCCCAGTTACCAGGTTAACCT')),
            Adapter('Native barcoding 4',
                    both_ends_sequence=('NB04', 'GGTGCTGTTCGGATTCTATCGTGTTTCCCTATTAACCT')),
            Adapter('Native barcoding 5',
                    both_ends_sequence=('NB05', 'GGTGCTGCTTGTCCAGGGTTTGTGTAACCTTTTAACCT')),
            Adapter('Native barcoding 6',
                    both_ends_sequence=('NB06', 'GGTGCTGTTCTCGCAAAGGCAGAAAGTAGTCTTAACCT')),
            Adapter('Native barcoding 7',
                    both_ends_sequence=('NB07', 'GGTGCTGGTGTTACCGTGGGAATGAATCCTTTTAACCT')),
            Adapter('Native barcoding 8',
                    both_ends_sequence=('NB08', 'GGTGCTGTTCAGGGAACAAACCAAGTTACGTTTAACCT')),
            Adapter('Native barcoding 9',
                    both_ends_sequence=('NB09', 'GGTGCTGAACTAGGCACAGCGAGTCTTGGTTTTAACCT')),
            Adapter('Native barcoding 10',
                    both_ends_sequence=('NB10', 'GGTGCTGAAGCGTTGAAACCTTTGTCCTCTCTTAACCT')),
            Adapter('Native barcoding 11',
                    both_ends_sequence=('NB11', 'GGTGCTGGTTTCATCTATCGGAGGGAATGGATTAACCT')),
            Adapter('Native barcoding 12',
                    both_ends_sequence=('NB12', 'GGTGCTGCAGGTAGAAAGAAGCAGAATCGGATTAACCT')),
            Adapter('PCR barcoding 1',
                    both_ends_sequence=('BC01', 'AAGAAAGTTGTCGGTGTCTTTGTG')),
            Adapter('PCR barcoding 2',
                    both_ends_sequence=('BC02', 'TCGATTCCGTTTGTAGTCGTCTGT')),
            Adapter('PCR barcoding 3',
                    both_ends_sequence=('BC03', 'GAGTCTTGTGTCCCAGTTACCAGG')),
            Adapter('PCR barcoding 4',
                    both_ends_sequence=('BC04', 'TTCGGATTCTATCGTGTTTCCCTA')),
            Adapter('PCR barcoding 5',
                    both_ends_sequence=('BC05', 'CTTGTCCAGGGTTTGTGTAACCTT')),
            Adapter('PCR barcoding 6',
                    both_ends_sequence=('BC06', 'TTCTCGCAAAGGCAGAAAGTAGTC')),
            Adapter('PCR barcoding 7',
                    both_ends_sequence=('BC07', 'GTGTTACCGTGGGAATGAATCCTT')),
            Adapter('PCR barcoding 8',
                    both_ends_sequence=('BC08', 'TTCAGGGAACAAACCAAGTTACGT')),
            Adapter('PCR barcoding 9',
                    both_ends_sequence=('BC09', 'AACTAGGCACAGCGAGTCTTGGTT')),
            Adapter('PCR barcoding 10',
                    both_ends_sequence=('BC10', 'AAGCGTTGAAACCTTTGTCCTCTC')),
            Adapter('PCR barcoding 11',
                    both_ends_sequence=('BC11', 'GTTTCATCTATCGGAGGGAATGGA')),
            Adapter('PCR barcoding 12',
                    both_ends_sequence=('BC12', 'CAGGTAGAAAGAAGCAGAATCGGA'))
            ]
