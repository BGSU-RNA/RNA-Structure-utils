import unittest

from rnastructure.tertiary.cif import CIF
from rnastructure.tertiary.cif import MissingColumn
from rnastructure.tertiary.cif import MissingBlockException

with open('test/files/1FAT.cif', 'rb') as raw:
    DATA = CIF(raw)


class SimpleCIFTest(unittest.TestCase):
    def setUp(self):
        self.cif = DATA

    def test_gets_table_with_leading_underscore(self):
        ans = 'pdbx_poly_seq_scheme'
        val = self.cif.table('_' + ans)
        self.assertEqual(val.name, ans)

    def test_gets_table_without_leading_underscore(self):
        ans = 'pdbx_poly_seq_scheme'
        val = self.cif.table(ans)
        self.assertEqual(val.name, ans)

    def test_attribute_gives_table(self):
        ans = 'pdbx_poly_seq_scheme'
        val = self.cif.pdbx_poly_seq_scheme
        self.assertEqual(val.name, ans)

    def test_fails_getting_unknown_table(self):
        with self.assertRaises(MissingBlockException):
            self.cif.table('bob')

    def test_raises_key_when_getting_missing_attribute(self):
        with self.assertRaises(AttributeError):
            self.cif.bob


class SimpleTableTest(unittest.TestCase):
    def setUp(self):
        self.data = DATA.table('pdbx_poly_seq_scheme')

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
        val = self.data.row(0)
        self.assertEqual(val, ans)

    def test_fails_getting_too_large_row(self):
        with self.assertRaises(IndexError):
            self.data.row(9000)

    def test_iterates_over_all_rows(self):
        ans = 1008
        val = len(list(self.data.rows))
        self.assertEqual(val, ans)

    def test_gets_a_column(self):
        ans = ['A', 'B', 'C', 'D']
        val = sorted(list(set(self.data.column('asym_id'))))
        self.assertEqual(val, ans)

    def test_fails_getting_missing_column(self):
        with self.assertRaises(MissingColumn):
            self.data.column('bob')

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

    def test_get_item_can_give_range(self):
        ans = [{'asym_id': 'A',
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
                'hetero': 'n'},
               {'asym_id': 'A',
                'entity_id': '1',
                'seq_id': '2',
                'mon_id': 'ASN',
                'ndb_seq_num': '2',
                'pdb_seq_num': '2',
                'auth_seq_num': '2',
                'pdb_mon_id': 'ASN',
                'auth_mon_id': 'ASN',
                'pdb_strand_id': 'A',
                'pdb_ins_code': '.',
                'hetero': 'n'}]
        val = self.data[0:2]
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
        with self.assertRaises(KeyError):
            self.data['bob']

    def test_dot_on_missing_column_gives_attribute(self):
        with self.assertRaises(AttributeError):
            self.data.bob

    def test_get_item_on_too_big_int_gives_index(self):
        with self.assertRaises(IndexError):
            self.data[90000]
