from __future__ import with_statement

import unittest

from rnastructure.tertiary.cif import CIF
# from rnastructure.tertiary.cif import Chain
from rnastructure.tertiary.cif import MissingColumn
from rnastructure.tertiary.cif import MissingBlockException


class SimpleCIFTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open('files/1FAT.cif', 'rb') as raw:
            cls.data = CIF(raw)

    def setUp(self):
        self.cif = self.__class__.data

    def test_gets_table_with_leading_underscore(self):
        val = self.cif.table('_pdbx_poly_seq_scheme')
        self.assertTrue(val is not None)

    def test_gets_table_without_leading_underscore(self):
        val = self.cif.table('pdbx_poly_seq_scheme')
        self.assertTrue(val is not None)

    def test_attribute_gives_table(self):
        val = self.cif.pdbx_poly_seq_scheme
        self.assertTrue(val is not None)

    def test_fails_getting_unknown_table(self):
        self.assertRaises(MissingBlockException, self.cif.table, 'bob')

    def test_raises_key_when_getting_missing_attribute(self):
        self.assertRaises(AttributeError, lambda: self.cif.bob)


class SimpleTableTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open('files/1FAT.cif', 'rb') as raw:
            cls.data = CIF(raw)

    def setUp(self):
        self.data = self.__class__.data.table('pdbx_poly_seq_scheme')

    def test_gets_all_columns(self):
        val = self.data.columns
        ans = ['asym_id', 'entity_id', 'seq_id', 'mon_id', 'ndb_seq_num',
               'pdb_seq_num', 'auth_seq_num', 'pdb_mon_id', 'auth_mon_id',
               'pdb_strand_id', 'pdb_ins_code', 'hetero']
        self.assertEquals(val, ans)

    def test_len_is_row_count(self):
        val = len(self.data)
        ans = 1008
        self.assertEqual(val, ans)

    def test_gets_size(self):
        ans = (1008, 12)
        val = self.data.size()
        self.assertEqual(val, ans)

    def test_gets_a_row(self):
        ans = {
            'asym_id': 'A',
            'entity_id': '1',
            'seq_id': '1',
            'mon_id': 'SER',
            'ndb_seq_num': '1',
            'pdb_seq_num': '1',
            'auth_seq_num': '1',
            'pdb_mon_id': 'SER',
            'auth_mon_id': 'SER',
            'pdb_strand_id': 'A',
            'pdb_ins_code': '.',
            'hetero': 'n'
        }
        val = self.data.rows[0]
        self.assertEqual(val, ans)

    def test_fails_getting_too_large_row(self):
        self.assertRaises(IndexError, lambda: self.data.rows[9000])

    def test_iterates_over_all_rows(self):
        ans = 1008
        val = len(self.data.rows)
        self.assertEqual(val, ans)

    def test_gets_a_column(self):
        ans = ['A', 'B', 'C', 'D']
        val = sorted(list(set(self.data.column('asym_id'))))
        self.assertEqual(val, ans)

    def test_fails_getting_missing_column(self):
        self.assertRaises(MissingColumn, self.data.column, 'bob')

    def test_get_item_can_give_row(self):
        ans = {
            'asym_id': 'A',
            'entity_id': '1',
            'seq_id': '1',
            'mon_id': 'SER',
            'ndb_seq_num': '1',
            'pdb_seq_num': '1',
            'auth_seq_num': '1',
            'pdb_mon_id': 'SER',
            'auth_mon_id': 'SER',
            'pdb_strand_id': 'A',
            'pdb_ins_code': '.',
            'hetero': 'n'
        }
        val = self.data[0]
        self.assertEqual(val, ans)

    def test_get_item_can_give_subtable(self):
        ans = ['SER', 'ASN']
        val = self.data[0:2].mon_id
        self.assertEqual(val, ans)

    def test_get_item_can_give_column(self):
        ans = ['A', 'B', 'C', 'D']
        val = sorted(list(set(self.data['asym_id'])))
        self.assertEqual(val, ans)

    def test_dot_gives_column(self):
        ans = ['A', 'B', 'C', 'D']
        val = sorted(list(set(self.data.asym_id)))
        self.assertEqual(val, ans)

    def test_get_item_on_missing_string_gives_key(self):
        self.assertRaises(KeyError, lambda: self.data['bob'])

    def test_dot_on_missing_column_gives_attribute(self):
        self.assertRaises(AttributeError, lambda: self.data.bob)

    def test_get_item_on_too_big_int_gives_index(self):
        self.assertRaises(IndexError, lambda: self.data[90000])


class SimpleSymmetryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open('files/1FAT.cif', 'rb') as raw:
            cls.data = CIF(raw)

    def setUp(self):
        self.cif = self.__class__.data
        self.operators = list(self.cif.symmetry_operators())

    def test_it_gets_all_operators(self):
        val = len(self.operators)
        ans = 1
        self.assertEquals(val, ans)

    def test_it_can_map_asym_to_symmetry(self):
        val = [op['name'] for op in self.cif.operators('A')]
        ans = ['1_555']
        self.assertEqual(val, ans)

    def test_it_gets_one_model(self):
        val = self.operators[0].model(1)
        self.assertTrue(val)

    def test_it_gets_model_with_string_number(self):
        val = self.operators[0].model('1')
        self.assertTrue(val)

    def test_it_gives_none_for_missing_model(self):
        val = self.operators[0].model('bob')
        self.assertTrue(val is None)

    def test_it_gets_all_models(self):
        val = len(list(self.operators[0].models()))
        ans = 1
        self.assertEqual(ans, val)


class SimpleModelTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open('files/1FAT.cif', 'rb') as raw:
            cls.data = CIF(raw)

    def setUp(self):
        self.cif = self.__class__.data
        self.models = list(self.cif.models())

    def test_it_finds_all_models(self):
        val = len(self.models)
        ans = 1
        self.assertEqual(ans, val)

    def test_it_can_get_one_chain(self):
        chain = self.models[0].chain('A')
        self.assertTrue(chain)

    def test_it_gives_none_for_missing_chain(self):
        chain = self.models[0].chain('bob')
        self.assertFalse(chain)

    def test_it_can_get_all_chains(self):
        chains = list(self.models[0].chains())
        val = len(chains)
        ans = 4
        self.assertEqual(ans, val)


class SimpleChainTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open('files/1FAT.cif', 'rb') as raw:
            cls.data = CIF(raw)

    def setUp(self):
        self.cif = self.__class__.data
        self.chains = list(self.cif.chains())

    def test_knows_if_it_is_true(self):
        self.assertTrue(self.chains[0])

    def test_can_get_all_residues(self):
        val = len(list(self.chains[0].residues()))
        ans = 232
        self.assertEqual(ans, val)

    def tes_it_knows_if_it_has_breaks(self):
        self.assertTrue(self.chains[0].has_breaks())


class SimplePolymerChainTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open('files/1FAT.cif', 'rb') as raw:
            cls.data = CIF(raw)

    def setUp(self):
        self.cif = self.__class__.data
        self.chain = self.cif.chain('1_555', '1', 'A')

    def test_unobs_is_sorted(self):
        val = sorted(self.chain.unobs, key=lambda u: int(u['auth_seq_id']))
        ans = [u for u in self.chain.unobs]
        self.assertEqual(ans, val)

    def test_atoms_are_sorted(self):
        val = sorted(self.chain.atoms(), key=lambda a: int(a['label_seq_id']))
        ans = [atom for atom in self.chain.atoms()]
        self.assertEqual(ans, val)

    def test_can_get_polymers(self):
        val = [len(p) for p in self.chain.polymers()]
        ans = [36, 232 - 36]
        self.assertEqual(ans, val)

    def test_it_has_no_breaks(self):
        self.assertFalse(self.chain.polymer(0).has_breaks())
        self.assertFalse(self.chain.polymer(1).has_breaks())


class CifPolymersWithInsCodeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open('files/2UUA.cif', 'rb') as raw:
            cls.data = CIF(raw)

    def setUp(self):
        self.cif = self.__class__.data
        self.data = list(self.cif.chain('1_555', '1', 'A').polymers())

    def test_can_get_chain_breaks_with_ins_code(self):
        val = [poly["chain"] for poly in self.data]
        ans = ['A', 'A']
        for poly in self.data:
            print(poly.residue(0))
            print(poly.residue(-1))
            if len(poly) < 10:
                print(poly.sequence)
            print(len(poly))
        self.assertEqual(ans, val)


class CifPolymersTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open('files/1FAT.cif', 'rb') as raw:
            cls.data = CIF(raw)

    def setUp(self):
        self.cif = self.__class__.data
        self.data = self.cif.chain('1_555', '1', 'D').polymers()

    def test_gets_all_polymers_in_all_chains(self):
        val = sorted(set([poly['chain'] for poly in self.cif.polymers()]))
        ans = ['A', 'B', 'C', 'D']
        self.assertEqual(val, ans)

    def test_gets_requested_chain_polymer(self):
        val = [poly['chain'] for poly in self.data]
        ans = ['D', 'D']
        self.assertEqual(val, ans)

    def test_finds_polymer_sequence(self):
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


# class SimpleResidueTest(unittest.TestCase):
#     @classmethod
#     def setUpClass(cls):
#         with open('files/1FAT.cif', 'rb') as raw:
#             cls.data = CIF(raw)

#     def setUp(self):
#         self.cif = self.__class__.data
#         self.residues = list(self.cif.residues())

    # def test_can_generate_a_unit_id(self):
    #     self.fail()
