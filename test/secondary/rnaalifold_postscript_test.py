import unittest

from rnastructure.secondary.rnaalifold_postscript import Parser


class RNAalifoldPostscriptTest(unittest.TestCase):
    def setUp(self):
        with open('test/files/alirna.ps', 'r') as raw:
            self.parser = Parser(raw)

    def test_load_sequence(self):
        val = self.parser.sequence
        ans = "GGGCCCGUAGCACAGUGGA__AGAGCACAUGCCUUCCAACCA_" + \
              "GAUGUCCCGGGUUCGAAUCCAGCCGAGCCCA"
        self.assertEquals(val, ans)
