import os
import unittest
from StringIO import StringIO

from rnastructure.primary.fold import Mfold
from rnastructure.primary.fold import UNAFold
from rnastructure.primary.fold import RNAalifold
from rnastructure.util.wrapper import InvalidInputError
from rnastructure.util.wrapper import ProgramTimeOutError


class RNAalifoldGenerateFileTest(unittest.TestCase):
    def setUp(self):
        self.fold = RNAalifold()
        self.sequences = ["ggggggggggggaaaaaaaacccccccccccc", 
                          "ggggggggggggaaaaaaaacccccccccccc"]
        io = StringIO()
        self.fold._generate_input_file_(io, self.sequences)
        self.contents = io.getvalue().split("\n")
        io.close()

    def test_header(self):
        val = self.contents[0]
        self.assertEqual(val, '# STOCKHOLM 1.0')

    def test_file_len(self):
        val = len(self.contents)
        self.assertEqual(val, 4)


class RNAalifoldTest(unittest.TestCase):
    def setUp(self):
        self.fold = RNAalifold()
        self.sequences = ["ggggggggggggaaaaaaaacccccccccccc", 
                          "ggggggggggggaaaaaaaacccccccccccc"]
        self.results = self.fold(self.sequences)

    def test_program(self):
        self.assertEqual(self.fold.program, 'RNAalifold')

    def test_result_count(self):
        val = len(self.results)
        self.assertEqual(val, 1)

    def test_result_sequence(self):
        val = self.results[0].sequence
        ans = 'GGGGGGGGGGGGAAAAAAAACCCCCCCCCCCC'
        self.assertEqual(val, ans)

    def test_result_indices(self):
        val = self.results[0].indices()
        ans = {'hairpin': [([12, 13, 14, 15, 16, 17, 18, 19],)]}
        self.assertEqual(val, ans)

    def test_result_loops(self):
        val = self.results[0].loops()
        ans = {'hairpin': ['AAAAAAAA']}
        self.assertEqual(val, ans)


class BasicUNAFoldTest(unittest.TestCase):

    def setUp(self):
        self.fold = UNAFold()
        self.sequence = "ggggggggggggaaaaaaaacccccccccccc"
        self.results = self.fold(self.sequence)

    def test_program(self):
        self.assertEqual(self.fold.program, 'UNAFold.pl')

    def test_max_length(self):
        folder = UNAFold(length=10)
        seq = 'a' * 11
        self.assertRaises(InvalidInputError, folder, seq)

    def test_timeout(self):
        folder = UNAFold(time=0.1, length=100)
        seq = 'cccccccccccccccccccccccaaaaaaaaaaaaaggggggggggggggggggggg'
        self.assertRaises(ProgramTimeOutError, folder, seq)

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


class BasicMfoldTest(unittest.TestCase):

    def setUp(self):
        self.fold = Mfold()
        self.sequence = "ggggggggggggaaaaaaaacccccccccccc"
        self.results = self.fold(self.sequence)

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
        fold = UNAFold()
        self.sequence = "ggggggggggggaaaaaaaacccccccccccc"
        self.result = fold(self.sequence)[0]

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
