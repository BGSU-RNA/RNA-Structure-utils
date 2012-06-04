import unittest

from rnastructure.secondary.pseudoknot import RemovePseudoknots
from rnastructure.secondary.dot_bracket import Parser as DotParser


class SimpleRemoveTest(unittest.TestCase):
    def setUp(self):
        self.remover = RemovePseudoknots()
        self.structure = '((..{{..))}}'
        self.raw = DotParser(self.structure)
        self.raw.sequence = 'ggaaccttccgg'
        self.parsed = self.remover(self.raw)

    def test_simple_removal(self):
        ans_parser = DotParser('((......))..')
        ans = ans_parser._pairs
        val = self.parsed._pairs
        self.assertEqual(val, ans)
