import unittest


import rnastructure.util.unit_ids as uid


class BaseGeneratorTest(unittest.TestCase):

    base = {'pdb': '1GID', 'type': 'BA1', 'model': 1, 'chain': 'A',
            'residue': 'C', 'number': 10}

    def simple_id(self):
        return self.generator(self.base)

    def non_standard_symmetry(self):
        data = dict(self.base, symmetry_operator='2_555')
        return self.generator(data)

    def insertion_id(self):
        data = dict(self.base, insertion_code='g')
        return self.generator(data)

    def qmark_insertion(self):
        data = dict(self.base)
        data['insertion_code'] = '?'
        return self.generator(data)

    def missing(self, *keys, **kwargs):
        data = dict(self.base)
        for key in keys:
            data.pop(key)
        data.update(kwargs)
        return self.generator(data)


class MatlabIdGeneratorTest(BaseGeneratorTest):

    def setUp(self):
        self.generator = uid.MatlabIdGenerator()

    def test_generates_simple_id(self):
        self.assertEqual(self.simple_id(), 'A:C10')

    def test_generates_with_insertion(self):
        self.assertEqual(self.insertion_id(), 'A:C10g')

    def test_treats_Qmark_insertion_as_empty(self):
        self.assertEqual(self.qmark_insertion(), 'A:C10')

    def test_fails_without_chain(self):
        with self.assertRaises(uid.MissingRequiredFragment):
            self.missing('chain')

    def test_fails_without_residue(self):
        with self.assertRaises(uid.MissingRequiredFragment):
            self.missing('residue')

    def test_fails_without_number(self):
        with self.assertRaises(uid.MissingRequiredFragment):
            self.missing('number')

    def test_succeeds_without_pdb(self):
        self.assertEqual(self.missing('pdb'), 'A:C10')


class NucleotideIdGeneratorTest(BaseGeneratorTest):

    def setUp(self):
        self.generator = uid.NucleotideIdGenerator()

    def test_generates_simple_id(self):
        self.assertEquals(self.simple_id(), '1GID_BA1_1_A_10_C_')

    def test_generates_with_insertion(self):
        self.assertEquals(self.insertion_id(), '1GID_BA1_1_A_10_C_g')

    def test_generates_with_Qmark_insertion_as_empty(self):
        self.assertEquals(self.qmark_insertion(), '1GID_BA1_1_A_10_C_')

    def test_fails_without_pdb(self):
        with self.assertRaises(uid.MissingRequiredFragment):
            self.missing('pdb')

    def test_fails_without_model(self):
        with self.assertRaises(uid.MissingRequiredFragment):
            self.missing('model')

    def test_fails_without_type(self):
        with self.assertRaises(uid.MissingRequiredFragment):
            self.missing('type')

    def test_fails_without_chain(self):
        with self.assertRaises(uid.MissingRequiredFragment):
            self.missing('chain')

    def test_fails_without_residue(self):
        with self.assertRaises(uid.MissingRequiredFragment):
            self.missing('residue')

    def test_fails_without_number(self):
        with self.assertRaises(uid.MissingRequiredFragment):
            self.missing('number')


class UnitIdGeneratorTest(BaseGeneratorTest):
    def setUp(self):
        self.generator = uid.UnitIdGenerator()

    def test_generates_simple_id(self):
        self.assertEquals(self.simple_id(), '1GID|1|A|C|10')

    def test_generates_with_insertion(self):
        self.assertEquals(self.insertion_id(), '1GID|1|A|C|10|||g')

    def test_generates_with_Qmark_insertion_as_empty(self):
        self.assertEquals(self.qmark_insertion(), '1GID|1|A|C|10')

    def test_generates_with_different_symmetry(self):
        self.assertEquals(self.non_standard_symmetry(),
                          '1GID|1|A|C|10||||2_555')

    def test_fails_without_pdb(self):
        with self.assertRaises(uid.MissingRequiredFragment):
            self.missing('pdb')

    def test_fails_without_model(self):
        with self.assertRaises(uid.MissingRequiredFragment):
            self.missing('model')

    def test_fails_without_chain(self):
        with self.assertRaises(uid.MissingRequiredFragment):
            self.missing('chain')

    def test_fails_without_number(self):
        with self.assertRaises(uid.MissingRequiredFragment):
            self.missing('number')

    def test_fails_without_residue(self):
        with self.assertRaises(uid.MissingRequiredFragment):
            self.missing('residue')

    def test_works_without_residue_and_number(self):
        self.assertEquals(self.missing('residue', 'number'), '1GID|1|A')

    def test_works_without_residue_and_number_and_chain(self):
        self.assertEquals(self.missing('residue', 'number', 'chain'), '1GID|1')

    def test_works_without_residue_and_number_and_chain_model(self):
        self.assertEquals(self.missing('residue', 'number', 'chain', 'model'),
                          '1GID')

    def test_fails_with_missing_many_and_non_standard_symmetry(self):
        with self.assertRaises(uid.MissingRequiredFragment):
            self.missing('residue', 'number', 'chain', 'model',
                         symmetry_operator='2_555')
