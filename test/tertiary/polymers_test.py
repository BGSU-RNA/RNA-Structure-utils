import unittest
from pprint import pprint

from Bio.PDB.Polypeptide import three_to_one

from rnastructure.tertiary.cif import CIF


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
        ans = [35, 196]
        self.assertEqual(val, ans)

    def test_finds_polymers(self):
        ans = [
            'SNDIYFNFQRFNETNLILQRDASVSSSGQLRLTNL',
            'NGEPRVGSLGRAFYSAPIQIWDNTTGTVASFATSFTFNIQVPNNAGPADGLAFALVPVGSQPK' +
            'DKGGFLGLFDGSNSNFHTVAVEFDTLYNKDWDPTERHIGIDVNSIRSIKTTRWDFVNGENAEV' +
            'LITYDSSTNLLVASLVYPSQKTSFIVSDTVDLKSVLPEWVSVGFSATTGINKGNVETNDVLSW' +
            'SFASKLS'
        ]
        val = []
        for poly in self.data:
            val.append(''.join([three_to_one(s) for s in poly.sequence]))
        pprint(val[1])
        pprint(ans[1])
        self.assertEqual(val, ans)
