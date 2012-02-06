import unittest
from rnastructure.secondary.basic import Node


class NodeTest(unittest.TestCase):
    def setUp(self):
        self.root = Node((None, 21))

    def test_eq(self):
      val = Node((None, 21))
      self.assertTrue(val == self.root)

    def test_add(self):
        child = Node(1, 2)
        self.root.add_child(child)
        val = self.root.children
        ans = [child]
        self.assertEqual(val, ans)

    def test_get_missing_largest(self):
        ans = Node((None, None))
        val = self.root.largest()
        self.assertEqual(val, ans)

    def test_spans(self):
        val = self.root.spans()
        ans = range(0, 21)
        self.assertEqual(val, ans)

    def test_unpaired_hairpin(self):
        val = self.root.unpaired()
        ans = tuple([range(0, 21)])
        self.assertEqual(val, ans)

    def test_nested_junction_unpaired(self):
        left = Node((3, 12))
        left.add_child(Node((5, 10)))
        right = Node((15, 20))
        right.add_child(Node((17, 19)))
        self.root.add_child(left)
        self.root.add_child(right)
        val = self.root.unpaired()
        ans = ([0, 1, 2], [13, 14])
        self.assertEqual(val, ans)

    def test_unpaired_junction(self):
        self.root.add_child(Node((4, 10)))
        self.root.add_child(Node((12, 18)))
        val = self.root.unpaired()
        ans = ([0, 1, 2, 3], [11], [19, 20])
        self.assertEqual(val, ans)

    def test_get_valid_largest(self):
        self.root.children.append(Node((1, 3)))
        ans = Node((1, 3))
        val = self.root.largest()
        self.assertEqual(val, ans)
