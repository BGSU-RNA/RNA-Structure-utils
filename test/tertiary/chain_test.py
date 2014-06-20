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

    def test_it_can_find_unit_id_for_index(self):
        val = self.data.residue(1).unit_id()
        ans = '1FAT|1|D|ASN|2'
        self.assertEqual(val, ans)
