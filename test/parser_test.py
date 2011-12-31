import unittest
from rnastructure import dot_bracket_parser as dbp


class OuterLoopTest(unittest.TestCase):
    def setUp(self):
        self.structure = "...((..((..))....)).."
        self.parser = dbp.Parser(self.structure)
        self.loops = self.parser.loops
        # self.sequences = self.parser.parse("...((..((..))....))..")
        # self.sequences = self.parser.parse("aaaccaaccaaggttttggtt")

    def test_hairpins(self):
        ans = [(9, 10)]
        self.assertEqual(self.loops['hairpins'], ans)

    def test_loops(self):
        ans = [([0, 1, 2], [19, 20]), ([5, 6], [13, 14, 15, 16])]
        self.assertEqual(self.loops['internal'], ans)


class SimpleTestDBPParser(unittest.TestCase):

    def setUp(self):
        self.structure = "((..((..))....))..((...))"
        self.parser = dbp.Parser(self.structure)
        self.loops = self.parser.loops

    def test_hairpins(self):
        ans = [(6, 7), (20, 21, 22)]
        self.assertEqual(self.loops['hairpins'], ans)

    def test_internal(self):
        ans = [([2, 3], [10, 11, 12, 13])]
        self.assertEqual(self.loops['internal'], ans)


class NestedTestDBPParser(unittest.TestCase):

    def setUp(self):
        self.structure = "...((..(((....)))...((((...))))...))..."
        self.parser = dbp.Parser(self.structure)
        self.loops = self.parser.loops

    def test_hairpins(self):
        ans = [(10, 11, 12, 13), (24, 25, 26)]
        self.assertEqual(self.loops['hairpins'], ans)
