import os
import unittest

from rnastructure.primary.fold import Mfold
from rnastructure.primary.fold import UNAfold


class BasicUNAfoldTest(unittest.TestCase):

    def setUp(self):
        self.fold = UNAfold()
        self.sequence = "ggggggggggggaaaaaaaacccccccccccc"
        self.results = self.fold.fold(self.sequence)

    def test_program(self):
        self.assertEqual(self.fold.program, 'UNAFold.pl')

    def test_result_count(self):
        val = len(self.results)
        self.assertEqual(val, 2)

    def test_get_result_sequence(self):
        first = self.results[0]
        val = first.sequence()
        ans = self.sequence
        self.assertEqual(val, ans)

    def test_get_result_pairing(self):
        parser = self.results[0].pairing()
        val = parser._pairs
        ans = [31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, None, None,
              None, None, None, None, None, None, 11, 10, 9, 8, 7, 6, 5, 4, 3,
              2, 1, 0]
        self.assertEqual(val, ans)


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
        val = first.sequence()
        ans = self.sequence
        self.assertEqual(val, ans)

    def test_get_result_pairing(self):
        parser = self.results[0].pairing()
        val = parser._pairs
        ans = [31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, None, None,
              None, None, None, None, None, None, 11, 10, 9, 8, 7, 6, 5, 4, 3,
              2, 1, 0]
        self.assertEqual(val, ans)


class ResultTest(unittest.TestCase):
    def setUp(self):
        fold = UNAfold()
        self.sequence = "ggggggggggggaaaaaaaacccccccccccc"
        self.result = fold.fold(self.sequence)[0]

    def test_indices(self):
        val = self.result.indices()
        ans = {'hairpin': [tuple([[12, 13, 14, 15, 16, 17, 18, 19]])]}
        self.assertEqual(val, ans)

    def test_loops(self):
        val = self.result.loops()
        ans = {'hairpin': ['aaaaaaaa']}
        self.assertEqual(val, ans)
