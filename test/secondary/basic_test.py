import unittest

from rnastructure.secondary.basic import Parser


class SimpleLoopTest(unittest.TestCase):
    def setUp(self):
        self.pairs = [None, None, 6, 5, None, 3, 2]
        self.parser = Parser(self.pairs)
        self.loops = self.parser.loops()
        self.flanking = self.parser.loops(flanking=True)

    def test_hairpins(self):
        ans = [tuple([[4]])]
        self.assertEqual(self.loops['hairpin'], ans)

    def test_flanking_hairpins(self):
        ans = [tuple([[3, 4, 5]])]
        self.assertEqual(self.flanking['hairpin'], ans)


class MultiLoopTest(unittest.TestCase):
    def setUp(self):
        self.pairs = [26, 25, None, None, 13, 12, None, 11, None, None, None,
                     7, 5, 4, None, 23, None, 21, None, None, None, 17, None,
                     15, None, 1, 0]
        self.parser = Parser(self.pairs)
        self.loops = self.parser.loops()

    def test_hairpin(self):
        val = self.loops['hairpin']
        ans = [tuple([[18, 19, 20]]), tuple([[8, 9, 10]])]
        self.assertEqual(val, ans)
