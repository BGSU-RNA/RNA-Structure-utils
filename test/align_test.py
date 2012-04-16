import unittest

from rnastructure.align import Align


class TrivalAlignTest(unittest.TestCase):
    def setUp(self):
        self.ref = "AAACC"
        self.align = Align(self.ref)
        self.result = self.align.align(self.ref)

    def test_ref_align(self):
        self.assertEqual(self.result['reference'], self.ref)

    def test_seq_align(self):
        self.assertEqual(self.result['sequence'], self.ref)


class SimpleAlignTest(unittest.TestCase):
    def setUp(self):
        self.ref = "AAACC"
        self.seq = "GGACC"
        self.align = Align(self.ref)
        self.result = self.align.align(self.seq)

    def test_ref_align(self):
        self.assertEqual(self.result['reference'], '--AAACC')

    def test_seq_align(self):
        self.assertEqual(self.result['sequence'], 'GGA--CC')
