import unittest

from rnastructure.secondary.rnaplot import Parser


class RNAalifoldPostscriptTest(unittest.TestCase):
    def setUp(self):
        with open('test/files/alirna.ps', 'r') as raw:
            self.parser = Parser(raw)

    def test_load_sequence(self):
        val = self.parser.ali_sequence
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
        val = self.parser.pairs
        ans = [72, 71, 70, 69, 68, 67, 66, None]
        self.assertEquals(val[0:8], ans)
