import unittest
from rnastructure.secondary.basic import Node


class NodeTest(unittest.TestCase):
    def setUp(self):
        self.root = Node((0, 20))

    def test_eq(self):
      val = Node((0, 20))
      self.assertTrue(val == self.root)

    def test_get_missing_largest(self):
        ans = Node((None, None))
        val = self.root.largest()
        self.assertEqual(val, ans)

    def test_get_valid_largest(self):
        self.root.children.append(Node((1, 3)))
        ans = Node((1, 3))
        val = self.root.largest()
        self.assertEqual(val, ans)

    # def test_add(self):
    #     child = Node(1, 2)
    #     self.root.add(child)
    #     val = self.root.children
    #     ans = [child]
    #     self.assertEqual(val, ans)
