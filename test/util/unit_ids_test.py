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


class GenericIdGeneratorTest(unittest.TestCase):
    def setUp(self):
        self.generator = uid.GenericIdGenerator({})

    def test_creates_none_insertion_code_for_qmark(self):
        val = self.generator.insertion_code({'insertion_code': '?'})
        self.assertTrue(val is None)

    def test_gives_none_insertion_code_for_spaces(self):
        val = self.generator.insertion_code({'insertion_code': ' '})
        self.assertTrue(val is None)

    def test_gives_none_insertion_code_for_empty(self):
        val = self.generator.insertion_code({'insertion_code': ''})
        self.assertTrue(val is None)

    def test_gives_none_for_missing_insertion_code(self):
        val = self.generator.insertion_code({})
        self.assertTrue(val is None)


class MatlabIdGeneratorTest(BaseGeneratorTest):

    def setUp(self):
        self.generator = uid.MatlabIdGenerator()

    def test_generates_simple_id(self):
        self.assertEqual(self.simple_id(), 'A:C10')

    def test_can_use_keyword_args(self):
        data = dict(self.base)
        data.pop('number')
        self.assertEquals(self.generator(data, number=3), 'A:C3')

    def test_given_takes_precedence_over_kwargs(self):
        data = dict(self.base)
        self.assertEquals(self.generator(data, number=3), 'A:C10')

    def test_generates_with_insertion(self):
        self.assertEqual(self.insertion_id(), 'A:C10g')

    def test_treats_Qmark_insertion_as_empty(self):
        self.assertEqual(self.qmark_insertion(), 'A:C10')

    def test_fails_without_chain(self):
        self.assertRaises(uid.MissingRequiredFragment, self.missing, 'chain')

    def test_fails_without_residue(self):
        self.assertRaises(uid.MissingRequiredFragment, self.missing, 'residue')

    def test_fails_without_number(self):
        self.assertRaises(uid.MissingRequiredFragment, self.missing, 'number')

    def test_succeeds_without_pdb(self):
        self.assertEqual(self.missing('pdb'), 'A:C10')


class NucleotideIdGeneratorTest(BaseGeneratorTest):

    def setUp(self):
        self.generator = uid.NucleotideIdGenerator()

    def test_generates_simple_id(self):
        self.assertEquals(self.simple_id(), '1GID_BA1_1_A_10_C_')

    def test_can_use_keyword_args(self):
        data = dict(self.base)
        data.pop('type')
        self.assertEquals(self.generator(data, type='AU'), '1GID_AU_1_A_10_C_')

    def test_given_takes_precedence_over_kwargs(self):
        self.assertEquals(self.generator(self.base, type='AU'),
                          '1GID_BA1_1_A_10_C_')

    def test_generates_with_insertion(self):
        self.assertEquals(self.insertion_id(), '1GID_BA1_1_A_10_C_g')

    def test_generates_with_Qmark_insertion_as_empty(self):
        self.assertEquals(self.qmark_insertion(), '1GID_BA1_1_A_10_C_')

    def test_fails_without_pdb(self):
        self.assertRaises(uid.MissingRequiredFragment, self.missing, 'pdb')

    def test_fails_without_model(self):
        self.assertRaises(uid.MissingRequiredFragment, self.missing, 'model')

    def test_fails_without_type(self):
        self.assertRaises(uid.MissingRequiredFragment, self.missing, 'type')

    def test_fails_without_chain(self):
        self.assertRaises(uid.MissingRequiredFragment, self.missing, 'chain')

    def test_fails_without_residue(self):
        self.assertRaises(uid.MissingRequiredFragment, self.missing, 'residue')

    def test_fails_without_number(self):
        self.assertRaises(uid.MissingRequiredFragment, self.missing, 'number')


class UnitIdGeneratorTest(BaseGeneratorTest):

    def setUp(self):
        self.generator = uid.UnitIdGenerator()

    def test_generates_simple_id(self):
        self.assertEquals(self.simple_id(), '1GID|1|A|C|10')

    def test_generates_using_kwargs(self):
        data = dict(self.base)
        data.pop('pdb')
        self.assertEquals(self.generator(data, pdb='2AW7'), '2AW7|1|A|C|10')

    def test_given_takes_precedence_over_kwargs(self):
        self.assertEquals(self.generator(self.base, pdb='2AW7'),
                          '1GID|1|A|C|10')

    def test_can_use_lookup_func(self):
        data = dict(self.base)
        data.pop('pdb')
        self.generator.lookup('pdb', lambda obj, **kw: '1J5E')
        val = self.generator(data)
        self.assertEqual(val, '1J5E|1|A|C|10')

    def test_generates_with_insertion(self):
        self.assertEquals(self.insertion_id(), '1GID|1|A|C|10|||g')

    def test_generates_with_Qmark_insertion_as_empty(self):
        self.assertEquals(self.qmark_insertion(), '1GID|1|A|C|10')

    def test_generates_with_different_symmetry(self):
        self.assertEquals(self.non_standard_symmetry(),
                          '1GID|1|A|C|10||||2_555')

    def test_fails_without_pdb(self):
        self.assertRaises(uid.MissingRequiredFragment, self.missing, 'pdb')

    def test_fails_without_model(self):
        self.assertRaises(uid.MissingRequiredFragment, self.missing, 'model')

    def test_fails_without_chain(self):
        self.assertRaises(uid.MissingRequiredFragment, self.missing, 'chain')

    def test_fails_without_number(self):
        self.assertRaises(uid.MissingRequiredFragment, self.missing, 'number')

    def test_fails_without_residue(self):
        self.assertRaises(uid.MissingRequiredFragment, self.missing, 'residue')

    def test_works_without_residue_and_number(self):
        self.assertEquals(self.missing('residue', 'number'), '1GID|1|A')

    def test_works_without_residue_and_number_and_chain(self):
        self.assertEquals(self.missing('residue', 'number', 'chain'), '1GID|1')

    def test_works_without_residue_and_number_and_chain_model(self):
        self.assertEquals(self.missing('residue', 'number', 'chain', 'model'),
                          '1GID')

    def test_fails_with_missing_many_and_non_standard_symmetry(self):
        self.assertRaises(uid.MissingRequiredFragment, self.missing, 'residue',
                          'number', 'chain', 'model',
                          symmetry_operator='2_555')


class GeneratedConverterTest(unittest.TestCase):
    def setUp(self):
        self.converter = uid.generate_converter('unit', 'nucleotide')

    def test_converts_using_arguments(self):
        val = self.converter('2AW7|1|A|G|11', type='AU')
        ans = '2AW7_AU_1_A_11_G_'
        self.assertEqual(ans, val)
