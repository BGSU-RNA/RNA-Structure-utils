import unittest

from rnastructure.secondary.dot_bracket import Parser
from rnastructure.secondary.dot_bracket import Writer


class RfamDialectTest(unittest.TestCase):
    def setUp(self):
        self.structure = '.AAA....<<<<aaa....>>>>'
        self.parser = Parser(self.structure, dialect='rfam')
        self.loops = self.parser.indices()

    def test_finds_hairpins(self):
        self.assertTrue('hairpin' in self.loops)

    def test_finds_no_internal(self):
        self.assertFalse('internal' in self.loops)


class OuterLoopTest(unittest.TestCase):
    def setUp(self):
        self.structure = "...((..((..))....)).."
        self.parser = Parser(self.structure)
        self.loops = self.parser.indices()
        self.flanking = self.parser.indices(flanking=True)

    def test_hairpins(self):
        ans = [tuple([[9, 10]])]
        self.assertEqual(self.loops['hairpin'], ans)

    def test_flanking_hairpins(self):
        ans = [tuple([[8, 9, 10, 11]])]
        self.assertEqual(self.flanking['hairpin'], ans)

    def test_external(self):
        ans = [([0, 1, 2], [19, 20])]
        self.assertEqual(self.loops['external'], ans)

    def test_flanking_external(self):
        ans = [([0, 1, 2, 3], [18, 19, 20])]
        self.assertEqual(self.flanking['external'], ans)

    def test_internal(self):
        ans = [([5, 6], [13, 14, 15, 16])]
        self.assertEqual(self.loops['internal'], ans)

    def test_flanking_internal(self):
        ans = [([4, 5, 6, 7], [12, 13, 14, 15, 16, 17])]
        self.assertEqual(self.flanking['internal'], ans)


class SimpleParserTest(unittest.TestCase):

    def setUp(self):
        self.structure = "((..((..))....))..((...))"
        self.parser = Parser(self.structure)
        self.loops = self.parser.indices()
        self.flanking = self.parser.indices(flanking=True)

    def test_hairpins(self):
        ans = [tuple([[6, 7]]), tuple([[20, 21, 22]])]
        self.assertEqual(self.loops['hairpin'], ans)

    def test_flanking_hairpins(self):
        ans = [tuple([[5, 6, 7, 8]]), tuple([[19, 20, 21, 22, 23]])]
        self.assertEqual(self.flanking['hairpin'], ans)

    def test_internal(self):
        ans = [([2, 3], [10, 11, 12, 13])]
        print(self.parser._tree.print_tree())
        self.assertEqual(self.loops['internal'], ans)

    def test_flanking_internal(self):
        ans = [([1, 2, 3, 4], [9, 10, 11, 12, 13, 14])]
        self.assertEqual(self.flanking['internal'], ans)


class NestedParserTest(unittest.TestCase):

    def setUp(self):
        self.structure = "...((..(((....)))...((((...))))...))..."
        self.parser = Parser(self.structure)
        self.loops = self.parser.indices()
        seq = "aaaccaacccuuuugggtttccccuuuggggtttggttt"
        self.sequences = self.parser.loops(seq)

    def test_hairpins(self):
        ans = [tuple([[24, 25, 26]]), tuple([[10, 11, 12, 13]])]
        self.assertEqual(self.loops['hairpin'], ans)

    def test_external_loop(self):
        ans = [([0, 1, 2], [36, 37, 38])]
        self.assertEqual(self.loops['external'], ans)

    def test_internal_loops(self):
        self.assertFalse('internal' in self.loops)

    # def test_junction_loops(self):
    #     ans = [([5, 6], [17, 18, 19], [31, 32, 33])]
    #     self.parser._tree.print_tree()
    #     self.assertEqual(self.loops['junction'], ans)

    def test_extract_hairpins(self):
        ans = ['uuu', 'uuuu']
        self.assertEqual(self.sequences['hairpin'], ans)

    # def test_extract_junction(self):
    #     ans = ['aa*ttt*ttt']
    #     self.assertEqual(self.sequences['junction'], ans)


# Commented out until I add pseudoknot parsing back in.
# class PseudoKnotTest(unittest.TestCase):
#     def setUp(self):
#         self.structure = ".((.{{..((..))..((..))...}}..))."
#         self.parser = Parser(self.structure)
#         self.loops = self.parser.indices()
#
#     def test_pseudoknots(self):
#         ans = [([4, 5], [25, 26])]
#         self.assertEqual(self.loops['pseudoknot'], ans)
#
#     def test_hairpins(self):
#         ans = [([10, 11]), ([18, 19])]
#         self.assertEqual(self.loops['hairpins'], ans)
#
#     def test_internal_loops(self):
#         ans = [([0], [31])]
#         self.assertEqual(self.loops['internal'], ans)
#
#     def test_junction_loops(self):
#         ans = [([3, 6, 7], [14, 15], [22, 23, 24, 27, 28])]
#         self.assertEqual(self.loops['junction'], ans)

class Ty1Test(unittest.TestCase):
    def setUp(self):
        self.structure = '(((((((.....((((((((...))))).)))......(((((((.........(((((.((.........(((((.(((((............(((.(((...((((....((..(((((....))))))).))))......((..(((...)))..)).))).)))........)))))))))).)).)))))..)))))))..)))).................(((.....))).((((....))))...((((((((.........(.(((((...((((..(((.........)))..))))....))))))))))))))........................((((((((....((((((..............))))))..............(((.(((((.((.......)).))))).)))....(((..(((((.((.((((........))))))((......))...........((((((((............))))))))......)))))..)))....(((((............))))).)))))..))).(((((......(((((((..........)))))))...)))))...................))).......'
        self.parser = Parser(self.structure)
        self.loops = self.parser.indices()

    def test_hairpins(self):
        ans = [tuple([[20,21,22]]),
               tuple([[121,122,123,124]]),
               tuple([[150,151,152]]),
               tuple([[230,231,232,233,234]]),
               tuple([[243,244,245,246]]),
               tuple([[290,291,292,293,294,295,296,297,298]]),
               tuple([[368,369,370,371,372,373,374,375,376,377,378,379,380,381]]),
               tuple([[414,415,416,417,418,419,420]]),
               tuple([[455,456,457,458,459,460,461,462]]),
               tuple([[471,472,473,474,475,476]]),
               tuple([[498,499,500,501,502,503,504,505,506,507,508,509]]),
               tuple([[543,544,545,546,547,548,549,550,551,552,553,554]]),
               tuple([[590,591,592,593,594,595,596,597,598,599]])]
        ans.reverse()
        self.assertEqual(self.loops['hairpin'], ans)


class SimpleWriterTest(unittest.TestCase):
    def setUp(self):
        self.structure = "...((..((..))....)).."
        self.parser = Parser(self.structure)
        self.formatter = Writer()

    def test_format(self):
        val = self.formatter.format(self.parser)
        self.assertEqual(val, self.structure)


class EmptyRightSideTest(unittest.TestCase):
    def setUp(self):
        self.structure = "((....((..))))"
        self.parser = Parser(self.structure)
        self.loops = self.parser.indices(flanking=True)

    def test_internal_loops(self):
        val = self.loops['internal']
        ans = [tuple([[1, 2, 3, 4, 5, 6], [11, 12]])]
        self.assertEqual(val, ans)


class EmptyLeftSideTest(unittest.TestCase):
    def setUp(self):
        self.structure = "((((..))....))"
        self.parser = Parser(self.structure)
        self.loops = self.parser.indices(flanking=True)

    def test_internal_loops(self):
        val = self.loops['internal']
        ans = [tuple([[1, 2], [7, 8, 9, 10, 11, 12]])]
        self.assertEqual(val, ans)
