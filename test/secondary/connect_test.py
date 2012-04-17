import unittest

from rnastructure.secondary.connect import Parser


class UnparserableConnectTest(unittest.TestCase):
    def test_complains(self):
        self.assertRaises(ValueError, Parser,
                          open('test/secondary/connect_test.py', 'r'))


class MultiStructureConnectTest(unittest.TestCase):
    def setUp(self):
        connect = open('test/files/simple_connect.ct', 'r')
        self.parser = Parser(connect)
        connect.close()
        self.loops = self.parser.loops()

    def test_parses_first_only(self):
        ans = [([9, 10, 11, 12, 13, 14, 15, 16],)]
        self.assertEqual(self.loops['hairpin'], ans)
