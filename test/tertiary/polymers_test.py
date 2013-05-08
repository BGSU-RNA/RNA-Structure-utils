import unittest

from rnastructure.tertiary.cif import CIF

from pprint import pprint


with open('test/files/1FAT.cif', 'rb') as raw:
    DATA = CIF(raw)


class PolymersTest(unittest.TestCase):
    def setUp(self):
        self.data = DATA.chain_polymer('D')

    def test_gets_all_chains(self):
        val = sorted(list(set([poly.chain for poly in DATA.polymers()])))
        ans = ['A', 'B', 'C', 'D']
        self.assertEqual(val, ans)

    def test_gets_polymers_with_breaks(self):
        val = [len(poly) for poly in self.data]
        ans = [34, 196]
        self.assertEqual(val, ans)
