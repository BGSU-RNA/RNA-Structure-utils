import unittest
from rnastructure import dot_bracket as dbp


class NodeTest(unittest.TestCase):

    def test_loading(self):
        node = dbp.Node()
        ans = 'a'
        node.add(ans)
        val = node[0]
        self.assertEqual(val, ans)

    def test_len(self):
        node = dbp.Node()
        node.add('c')
        node.add('.')
        self.assertEqual(len(node), 2)

    def test_iter(self):
        node = dbp.Node()
        node.add('a')
        node.add('b')
        node.add('c')
        ans = ['a', 'b', 'c']
        val = [c for c in node]
        self.assertEqual(val, ans)

    def test_empty_get_loop(self):
        node = dbp.Node()
        node.add(dbp.Node())
        ans = tuple()
        val = node.get_loop()
        self.assertEqual(val, ans)

    def test_simple_get_loop(self):
        node = dbp.Node()
        node.add(1)
        node.add(2)
        ans = tuple([[1, 2]])
        val = node.get_loop()
        self.assertEqual(val, ans)

    def test_internal_get_loop(self):
        node = dbp.Node()
        node.add(dbp.Node())
        node.add(1)
        node.add(2)
        node.add(dbp.Node())
        node.add(dbp.Node())
        node.add(3)
        ans = tuple([[1, 2], [3]])
        val = node.get_loop()
        self.assertEqual(val, ans)

    def test_outer_get_loop(self):
        node = dbp.Node()
        node.add(1)
        node.add(2)
        node.add(dbp.Node())
        node.add(3)
        ans = tuple([[1, 2], [3]])
        val = node.get_loop()
        self.assertEqual(val, ans)

    def test_junction_get_loop(self):
        node = dbp.Node()
        node.add(dbp.Node())
        node.add(1)
        node.add(2)
        node.add(dbp.Node())
        node.add(dbp.Node())
        node.add(3)
        node.add(4)
        node.add(5)
        node.add(dbp.Node())
        node.add(dbp.Node())
        node.add(6)
        node.add(7)
        ans = tuple([[1, 2], [3, 4, 5], [6, 7]])
        val = node.get_loop()
        self.assertEqual(val, ans)

