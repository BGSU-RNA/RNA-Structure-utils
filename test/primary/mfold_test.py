import os
import unittest

from rnastructure.primary.mfold import Mfold


class BasicMfoldTest(unittest.TestCase):

    def setUp(self):
        self.mfold = Mfold()
        self.sequence = "ggggggggggggaaaaaaaacccccccccccc"
        self.results = self.mfold.fold(self.sequence)

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
