import unittest

from rnastructure.tertiary.cif import CIF


class PolymersTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open('files/1FAT.cif', 'rb') as raw:
            cls.data = CIF(raw)

    def setUp(self):
        self.cif = self.__class__.data
        self.data = self.cif.chain('1_555', 1, 'D')

    def test_it_knows_the_chain_id(self):
        val = self.data['chain']
        ans = 'D'
        self.assertEqual(ans, val)

    def test_it_only_has_rows_from_correct_chain(self):
        val = list({a['auth_asym_id'] for a in self.data.atoms()})
        ans = ['D']
        self.assertEqual(ans, val)

    def test_it_knows_the_sequence(self):
        val = self.data.sequence
        ans = ['SER', 'ASN', 'ASP', 'ILE', 'TYR', 'PHE', 'ASN', 'PHE', 'GLN',
               'ARG', 'PHE', 'ASN', 'GLU', 'THR', 'ASN', 'LEU', 'ILE', 'LEU',
               'GLN', 'ARG', 'ASP', 'ALA', 'SER', 'VAL', 'SER', 'SER', 'SER',
               'GLY', 'GLN', 'LEU', 'ARG', 'LEU', 'THR', 'ASN', 'LEU', 'ASN',
               'GLY', 'ASN', 'GLY', 'GLU', 'PRO', 'ARG', 'VAL', 'GLY', 'SER',
               'LEU', 'GLY', 'ARG', 'ALA', 'PHE', 'TYR', 'SER', 'ALA', 'PRO',
               'ILE', 'GLN', 'ILE', 'TRP', 'ASP', 'ASN', 'THR', 'THR', 'GLY',
               'THR', 'VAL', 'ALA', 'SER', 'PHE', 'ALA', 'THR', 'SER', 'PHE',
               'THR', 'PHE', 'ASN', 'ILE', 'GLN', 'VAL', 'PRO', 'ASN', 'ASN',
               'ALA', 'GLY', 'PRO', 'ALA', 'ASP', 'GLY', 'LEU', 'ALA', 'PHE',
               'ALA', 'LEU', 'VAL', 'PRO', 'VAL', 'GLY', 'SER', 'GLN', 'PRO',
               'LYS', 'ASP', 'LYS', 'GLY', 'GLY', 'PHE', 'LEU', 'GLY', 'LEU',
               'PHE', 'ASP', 'GLY', 'SER', 'ASN', 'SER', 'ASN', 'PHE', 'HIS',
               'THR', 'VAL', 'ALA', 'VAL', 'GLU', 'PHE', 'ASP', 'THR', 'LEU',
               'TYR', 'ASN', 'LYS', 'ASP', 'TRP', 'ASP', 'PRO', 'THR', 'GLU',
               'ARG', 'HIS', 'ILE', 'GLY', 'ILE', 'ASP', 'VAL', 'ASN', 'SER',
               'ILE', 'ARG', 'SER', 'ILE', 'LYS', 'THR', 'THR', 'ARG', 'TRP',
               'ASP', 'PHE', 'VAL', 'ASN', 'GLY', 'GLU', 'ASN', 'ALA', 'GLU',
               'VAL', 'LEU', 'ILE', 'THR', 'TYR', 'ASP', 'SER', 'SER', 'THR',
               'ASN', 'LEU', 'LEU', 'VAL', 'ALA', 'SER', 'LEU', 'VAL', 'TYR',
               'PRO', 'SER', 'GLN', 'LYS', 'THR', 'SER', 'PHE', 'ILE', 'VAL',
               'SER', 'ASP', 'THR', 'VAL', 'ASP', 'LEU', 'LYS', 'SER', 'VAL',
               'LEU', 'PRO', 'GLU', 'TRP', 'VAL', 'SER', 'VAL', 'GLY', 'PHE',
               'SER', 'ALA', 'THR', 'THR', 'GLY', 'ILE', 'ASN', 'LYS', 'GLY',
               'ASN', 'VAL', 'GLU', 'THR', 'ASN', 'ASP', 'VAL', 'LEU', 'SER',
               'TRP', 'SER', 'PHE', 'ALA', 'SER', 'LYS', 'LEU', 'SER', 'ASP',
               'GLU', 'THR', 'THR', 'SER', 'GLU', 'GLY', 'LEU', 'ASN', 'LEU',
               'ALA', 'ASN', 'LEU', 'VAL', 'LEU', 'ASN', 'LYS', 'ILE', 'LEU']
        self.assertEqual(ans, val)

    def test_it_can_find_unit_id_for_index(self):
        val = self.data.unit_id(1)
        ans = '1FAT|1|D|ASN|1'
        self.assertEqual(val, ans)
