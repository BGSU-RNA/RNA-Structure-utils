import unittest

from StringIO import StringIO

from rnastructure.secondary import dot_bracket as DB

from rnastructure.secondary.connect import Parser
from rnastructure.secondary.connect import Writer
from rnastructure.secondary.connect import InvalidConnectLine


class UnparserableConnectTest(unittest.TestCase):
    def test_complains(self):
        self.assertRaises(InvalidConnectLine, Parser,
                          open('test/secondary/connect_test.py', 'r'))


class MultiStructureConnectTest(unittest.TestCase):
    def setUp(self):
        connect = open('test/files/simple_connect.ct', 'r')
        self.parser = Parser(connect)
        connect.close()
        self.loops = self.parser.indices()

    def test_parses_first_only(self):
        ans = [([9, 10, 11, 12, 13, 14, 15, 16],)]
        self.assertEqual(self.loops['hairpin'], ans)


class ConnectWriterTest(unittest.TestCase):
    def setUp(self):
        dot_parser = DB.Parser('((..))')
        dot_parser.sequence = 'ccaagg'
        self.writer = Writer()
        self.connect = self.writer.format(dot_parser)

    def test_formats_correctly(self):
        ans = ['6 Energy = ', 
               '1\tc\t0\t2\t6\t1',
               '2\tc\t1\t3\t5\t2',
               '3\ta\t2\t4\t0\t3',
               '4\ta\t3\t5\t0\t4',
               '5\tg\t4\t6\t2\t5',
               '6\tg\t5\t0\t1\t6',
               '']
        val = self.connect.split('\n')
        self.assertEqual(val, ans)

    def test_roundtrip(self):
        lines = StringIO(self.connect)
        parser = Parser(lines)
        val = parser.indices()['hairpin']
        ans = [([2, 3],)]
        self.assertEqual(val, ans)

class ConnectWriterTest(unittest.TestCase):
    def setUp(self):
        dot_parser = DB.Parser('((..))')
        self.writer = Writer()
        self.connect = self.writer.format(dot_parser)

    def test_formats_correctly(self):
        ans = ['6 Energy = ', 
               '1\t?\t0\t2\t6\t1',
               '2\t?\t1\t3\t5\t2',
               '3\t?\t2\t4\t0\t3',
               '4\t?\t3\t5\t0\t4',
               '5\t?\t4\t6\t2\t5',
               '6\t?\t5\t0\t1\t6',
               '']
        val = self.connect.split('\n')
        self.assertEqual(val, ans)

    def test_roundtrip(self):
        lines = StringIO(self.connect)
        parser = Parser(lines)
        val = parser.indices()['hairpin']
        ans = [([2, 3],)]
        self.assertEqual(val, ans)
