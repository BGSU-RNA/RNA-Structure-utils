from __future__ import with_statement

import unittest

# from rnastructure.secondary import dot_bracket
from rnastructure.secondary import rnaplot as rp


class PostScriptParserTest(unittest.TestCase):
    def setUp(self):
        with open('files/alirna.ps', 'r') as raw:
            self.parser = rp.PostScriptParser(raw)

    def test_load_sequence(self):
        val = self.parser.sequence
        ans = "GGGCCCGUAGCACAGUGGA__AGAGCACAUGCCUUCCAACCA_" + \
              "GAUGUCCCGGGUUCGAAUCCAGCCGAGCCCA"
        self.assertEquals(val, ans)

    def test_load_locations(self):
        val = self.parser.locations
        ans = [(110.62167358, 226.01321411),
               (109.99353790, 211.02636719),
               (109.36540222, 196.03953552)]
        self.assertEquals(len(val), 74)
        self.assertEquals(val[0:3], ans)

    def test_load_pairs(self):
        val = self.parser._pairs
        ans = [72, 71, 70, 69, 68, 67, 66, None]
        self.assertEqual(len(val), 74)
        self.assertEquals(val[0:8], ans)

    def test_load_bounding_box(self):
        val = self.parser.box
        ans = (66, 210, 518, 662)
        self.assertEqual(val, ans)


class SVGParserTest(unittest.TestCase):
    def setUp(self):
        with open('files/rna.svg', 'rb') as raw:
            self.parser = rp.SVGParser(raw)

    def test_loads_sequence(self):
        val = self.parser.sequence
        ans = 'AAAAAAAAAAAcccccccUUUUUUUUUUU'
        self.assertEqual(val, ans)

    def test_loads_locations(self):
        val = self.parser.locations
        ans = [(92.5, -101.743),
               (92.5, -86.743)]
        self.assertEqual(ans, val[0:2])

    def test_loads_pairs(self):
        val = self.parser._pairs
        ans = [28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, None, None, None,
               None, None, None, None, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
        self.assertEqual(ans, val)


# class WrapperTest(unittest.TestCase):
#     def setUp(self):
#         self.plotter = RNAplot()
#         self.data = dot_bracket.Parser("((((((((((......))))))))))")
#         self.data.sequence = "ccccccccccaaaaaagggggggggg"

#     def test_draws_using_svg(self):
#         val = self.plotter(self.data, output_format='svg')
#         self.fail()

#     def test_draws_using_ps(self):
#         val = self.plotter(self.data, output_format='ps')
#         self.fail()
