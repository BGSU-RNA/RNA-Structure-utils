from __future__ import with_statement

import unittest

from rnastructure.tertiary.cif import CIF


with open('files/1FAT.cif', 'rb') as raw:
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
        ans = [35, 196]
        self.assertEqual(val, ans)

    def test_finds_polymers(self):
        ans = [
            ['SER', 'ASN', 'ASP', 'ILE', 'TYR', 'PHE', 'ASN', 'PHE', 'GLN',
             'ARG', 'PHE', 'ASN', 'GLU', 'THR', 'ASN', 'LEU', 'ILE', 'LEU',
             'GLN', 'ARG', 'ASP', 'ALA', 'SER', 'VAL', 'SER', 'SER', 'SER',
             'GLY', 'GLN', 'LEU', 'ARG', 'LEU', 'THR', 'ASN', 'LEU'],
            ['ASN', 'GLY', 'GLU', 'PRO', 'ARG', 'VAL', 'GLY', 'SER', 'LEU',
             'GLY', 'ARG', 'ALA', 'PHE', 'TYR', 'SER', 'ALA', 'PRO', 'ILE',
             'GLN', 'ILE', 'TRP', 'ASP', 'ASN', 'THR', 'THR', 'GLY', 'THR',
             'VAL', 'ALA', 'SER', 'PHE', 'ALA', 'THR', 'SER', 'PHE', 'THR',
             'PHE', 'ASN', 'ILE', 'GLN', 'VAL', 'PRO', 'ASN', 'ASN', 'ALA',
             'GLY', 'PRO', 'ALA', 'ASP', 'GLY', 'LEU', 'ALA', 'PHE', 'ALA',
             'LEU', 'VAL', 'PRO', 'VAL', 'GLY', 'SER', 'GLN', 'PRO', 'LYS',
             'ASP', 'LYS', 'GLY', 'GLY', 'PHE', 'LEU', 'GLY', 'LEU', 'PHE',
             'ASP', 'GLY', 'SER', 'ASN', 'SER', 'ASN', 'PHE', 'HIS', 'THR',
             'VAL', 'ALA', 'VAL', 'GLU', 'PHE', 'ASP', 'THR', 'LEU', 'TYR',
             'ASN', 'LYS', 'ASP', 'TRP', 'ASP', 'PRO', 'THR', 'GLU', 'ARG',
             'HIS', 'ILE', 'GLY', 'ILE', 'ASP', 'VAL', 'ASN', 'SER', 'ILE',
             'ARG', 'SER', 'ILE', 'LYS', 'THR', 'THR', 'ARG', 'TRP', 'ASP',
             'PHE', 'VAL', 'ASN', 'GLY', 'GLU', 'ASN', 'ALA', 'GLU', 'VAL',
             'LEU', 'ILE', 'THR', 'TYR', 'ASP', 'SER', 'SER', 'THR', 'ASN',
             'LEU', 'LEU', 'VAL', 'ALA', 'SER', 'LEU', 'VAL', 'TYR', 'PRO',
             'SER', 'GLN', 'LYS', 'THR', 'SER', 'PHE', 'ILE', 'VAL', 'SER',
             'ASP', 'THR', 'VAL', 'ASP', 'LEU', 'LYS', 'SER', 'VAL', 'LEU',
             'PRO', 'GLU', 'TRP', 'VAL', 'SER', 'VAL', 'GLY', 'PHE', 'SER',
             'ALA', 'THR', 'THR', 'GLY', 'ILE', 'ASN', 'LYS', 'GLY', 'ASN',
             'VAL', 'GLU', 'THR', 'ASN', 'ASP', 'VAL', 'LEU', 'SER', 'TRP',
             'SER', 'PHE', 'ALA', 'SER', 'LYS', 'LEU', 'SER']
        ]
        val = [poly.sequence for poly in self.data]
        self.assertEqual(val, ans)
