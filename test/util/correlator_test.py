import unittest

from rnastructure.util import correlator as corr


class SimpleColumnCorrelationTest(unittest.TestCase):
    def setUp(self):
        self.seq = "-a--a-gg"
        self.correlations = corr.columns_to_index(self.seq)

    def test_correlations(self):
        ans = {1: 0, 4: 1, 6: 2, 7: 3}
        self.assertEqual(self.correlations, ans)


class NamedColumnCorrelationTest(unittest.TestCase):
    def setUp(self):
        self.seq = "aac---gg-t"
        self.names = ['1', '2', '3', '3A', '4', '5']
        self.correlations = corr.columns_to_index(self.seq, names=self.names)

    def test_correlations(self):
        ans = {0: '1', 1: '2', 2: '3', 6: '3A', 7: '4', 9: '5'}
        self.assertEqual(self.correlations, ans)


class ColumnCorrelationTest(unittest.TestCase):
    def setUp(self):
        self.seq = 'aaa--c-'
        self.ref = '-aa-ccu'
        self.correlations = corr.correlated_columns(self.ref, self.seq)

    def test_correlations(self):
        ans = [1, 2, 5]
        self.assertEqual(self.correlations, ans)


class NamedAlignedIndeciesCorrelationTest(unittest.TestCase):
    def setUp(self):
        seq = '-aa--cu'
        ref = 'uaag-c-'
        seq_names = ['a', 'b', 'c', 'd']
        ref_names = ['1', '1A', '2', '2A', '3']
        self.correlations = corr.correlate_aligned_indecies(ref, seq,
                                                            reference_names=ref_names,
                                                            sequence_names=seq_names)

    def test_correlations(self):
        ans = {'1A': 'a', '2': 'b', '3': 'c'}
        self.assertEqual(self.correlations, ans)
