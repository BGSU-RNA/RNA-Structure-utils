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
