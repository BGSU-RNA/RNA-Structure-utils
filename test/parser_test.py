import unittest
from rnastructure import dot_bracket as dbp


class OuterLoopTest(unittest.TestCase):
    def setUp(self):
        self.structure = "...((..((..))....)).."
        self.parser = dbp.Parser(self.structure)
        self.loops = self.parser.loops

    def test_hairpins(self):
        ans = [([9, 10])]
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
        ans = [([6, 7]), ([20, 21, 22])]
        self.assertEqual(self.loops['hairpins'], ans)

    def test_internal(self):
        ans = [([16, 17], []), ([2, 3], [10, 11, 12, 13])]
        self.assertEqual(self.loops['internal'], ans)


class NestedTestDBPParser(unittest.TestCase):

    def setUp(self):
        self.structure = "...((..(((....)))...((((...))))...))..."
        self.parser = dbp.Parser(self.structure)
        self.loops = self.parser.loops
        self.sequences = self.parser.parse("aaaccaacccuuuugggtttccccuuuggggtttggttt")

    def test_hairpins(self):
        ans = [([10, 11, 12, 13]), ([24, 25, 26])]
        self.assertEqual(self.loops['hairpins'], ans)

    def test_internal_loops(self):
        ans = [([0, 1, 2,], [36, 37, 38])]
        self.assertEqual(self.loops['internal'], ans)

    def test_extract_hairpins(self):
        ans = ['uuuu', 'uuu']
        self.assertEqual(self.sequences['hairpins'], ans)

    def test_extract_loops(self):
        ans = ['aaa*ttt']
        self.assertEqual(self.sequences['internal'], ans)
