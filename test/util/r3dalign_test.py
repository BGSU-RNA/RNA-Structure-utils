import unittest

import rnastructure.util.r3dalign as r3d


class R3DAlignBP2NTTest(unittest.TestCase):
    """The loader should read the csv file and produce a list of (n1, nt2)
    correlations from the csv file of interactions. We will attempt to make
    the output as close to unit ids as possible.
    """

    def setUp(self):
        with open('files/r3d.csv', 'rb') as raw:
            self.data = r3d.bp2nt(raw)

    def test_loads_all_data(self):
        self.assertEqual(3071, len(self.data))


class AsNtIdTest(unittest.TestCase):
    def test_handles_simple_id(self):
        ans = '1ABC|1|A|A|6'
        val = r3d.as_nt_id('1ABC', 'A:A6')
        self.assertEquals(ans, val)

    def test_handles_insertion_code(self):
        ans = '1ABC|1|A|A|2801|||A'
        val = r3d.as_nt_id('1ABC', 'A:A2801A')
        self.assertEquals(ans, val)
