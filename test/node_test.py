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
