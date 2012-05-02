import os
import unittest

from rnastructure.primary.fold import Mfold
from rnastructure.primary.fold import UNAfold
from rnastructure.primary.fold import FoldingFailedError
from rnastructure.primary.fold import FoldingTimeOutError


class BasicUNAfoldTest(unittest.TestCase):

    def setUp(self):
        self.fold = UNAfold()
        self.sequence = "ggggggggggggaaaaaaaacccccccccccc"
        self.results = self.fold.fold(self.sequence)

    def test_program(self):
        self.assertEqual(self.fold.program, 'UNAFold.pl')

    def test_max_length(self):
        folder = UNAfold(length=10)
        seq = 'a' * 11
        self.assertRaises(ValueError, folder.fold, seq)

    def test_timeout(self):
        folder = UNAfold(time=0.1, length=100)
        seq = 'cccccccccccccccccccccccaaaaaaaaaaaaaggggggggggggggggggggg'
        self.assertRaises(FoldingTimeOutError, folder.fold, seq)

    def test_result_count(self):
        val = len(self.results)
        self.assertEqual(val, 2)

    def test_get_result_sequence(self):
        first = self.results[0]
        val = first.sequence
        ans = self.sequence
        self.assertEqual(val, ans)

    def test_get_result_indices(self):
        val = self.results[0].parser._pairs
        ans = [31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, None, None,
              None, None, None, None, None, None, 11, 10, 9, 8, 7, 6, 5, 4, 3,
              2, 1, 0]
        self.assertEqual(val, ans)


# class LongFoldTest(unittest.TestCase):
#     def setUp(self):
#         self.sequence = 'GAGGAGAACTTCTAGTGTATATTCTGTATACCTAATATTATAGCCTTTATCAACAATGGAATCCCAACAATTATCTCAACATTCCCCGATT-TTTCATGGTAGCGCCTGTGCTTCGGTTACTTCTAAAGAAGTCCAAACAACTCAAGATCCGTTAGACATTTCAGCTTCCAAAACAGAAGAATGTGAGAAGGTTTCCACTCAGGCTAATTCTCAACAGCCAACAACACCTCCCTCATCTGCTGTTCCAGAGAACCATCATCATGCCTCTCCTCAAGCTGCTCAAGTACCATTG-CCACAAAATGGGCCGTACCCACAGCAGCGCATGATGAATACCCAA---CAAGCCAATATTTCTGGCTGGCCAGTATACGGGCACCC-ATCCTTGATGCCGTATCCACCTTATCAAATGTCACCTATGTACGCTCCACCTGGGGCACAATCACAGTTTACACAATATCCACAATATGTTGGAACACATTTGAACACCCCGTCACCTGAGTCAGGTAATTCATTTCCTGATTCATC-CTCAGCAAAGTCTAA---TATGACATCCACTAATCAACATGTCAGACCACCGCCAATCTTAACCTCACCTAATGACTTTCTAAATTGGGTTAAAATATACATCAAATTTTTACAAAATTCGAATCTC'
#         self.folder = UNAfold(length=len(self.sequence) + 1)
# 
#     def internal

class BasicMfoldTest(unittest.TestCase):

    def setUp(self):
        self.fold = Mfold()
        self.sequence = "ggggggggggggaaaaaaaacccccccccccc"
        self.results = self.fold.fold(self.sequence)

    def test_program(self):
        self.assertEqual(self.fold.program, 'mfold')

    def test_result_count(self):
        val = len(self.results)
        self.assertEqual(val, 2)

    def test_get_result_sequence(self):
        first = self.results[0]
        val = first.sequence
        ans = self.sequence
        self.assertEqual(val, ans)

    def test_get_result_indices(self):
        val = self.results[0].parser._pairs
        ans = [31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, None, None,
              None, None, None, None, None, None, 11, 10, 9, 8, 7, 6, 5, 4, 3,
              2, 1, 0]
        self.assertEqual(val, ans)


class ResultTest(unittest.TestCase):
    def setUp(self):
        fold = UNAfold()
        self.sequence = "ggggggggggggaaaaaaaacccccccccccc"
        self.result = fold.fold(self.sequence)[0]

    def test_get_connect_file(self):
        val = self.result.connect_file()
        ans = '32	dG = -31.8	sequence\n'
        self.assertEqual(val[0], ans)

    def test_indices(self):
        val = self.result.indices()
        ans = {'hairpin': [tuple([[12, 13, 14, 15, 16, 17, 18, 19]])]}
        self.assertEqual(val, ans)

    def test_loops(self):
        val = self.result.loops()
        ans = {'hairpin': ['aaaaaaaa']}
        self.assertEqual(val, ans)
