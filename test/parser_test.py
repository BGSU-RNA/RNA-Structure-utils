import unittest
from rna_structure_parser import RNAStructureParser


class TestParser(unittest.TestCase):
    """Test the RNA parser"""

    def loops(self, structure):
        self.parser = RNAStructureParser(structure)
        self.loops = self.parser.loops

    def test_simple_hairpins(self):
        self.loops("((..((..))..((..))..))..((..))")
        self.assertEqual(self.loops['hairpins'], [(2, 3)])
