import unittest

from rnastructure.secondary.basic import Parser
from rnastructure.secondary.basic import EmptyStructureError


class BadInputTest(unittest.TestCase):
    def test_empty_loops(self):
        self.assertRaises(EmptyStructureError, Parser, [])


class SimpleIndicesTest(unittest.TestCase):
    def setUp(self):
        self.pairs = [None, None, 6, 5, None, 3, 2]
        self.sequence = 'aacugcc'
        self.parser = Parser(self.pairs)
        self.indices = self.parser.indices()
        self.flanking = self.parser.indices(flanking=True)

    def test_hairpins(self):
        ans = [tuple([[4]])]
        self.assertEqual(self.indices['hairpin'], ans)

    def test_hairpin_sequences(self):
        ans = ['g']
        val = self.parser.loops(self.sequence)
        self.assertEqual(val['hairpin'], ans)

    def test_flanking_hairpins(self):
        ans = [tuple([[3, 4, 5]])]
        self.assertEqual(self.flanking['hairpin'], ans)

    def test_hairpin_flanking_sequences(self):
        ans = ['ugc']
        val = self.parser.loops(self.sequence, flanking=True)
        self.assertEqual(val['hairpin'], ans)


class MultiLoopTest(unittest.TestCase):
    def setUp(self):
        self.pairs = [26, 25, None, None, 13, 12, None, 11, None, None, None,
                     7, 5, 4, None, 23, None, 21, None, None, None, 17, None,
                     15, None, 1, 0]
        self.parser = Parser(self.pairs)
        self.indices = self.parser.indices()

    def test_hairpin(self):
        val = self.indices['hairpin']
        ans = [tuple([[18, 19, 20]]), tuple([[8, 9, 10]])]
        self.assertEqual(val, ans)
