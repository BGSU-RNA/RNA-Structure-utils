from Bio import pairwise2 as nw

from rnastructure.util import correlator as cor


class Align(object):
    """Class to align RNA sequences. """

    def __init__(self, reference=None):
        self.reference = reference
        self.match = 2
        self.mismatch = -1
        self.open = -0.5
        self.extend = -0.1

    def align_all(sequences):
        map(self.align, sequences)

    def best_alignment(self, sequences):
        alignments = self.align_all(sequences)
        return max(alignments, key=lambda a: a['score'])

    def align(self, sequence):
        alignments = nw.align.globalms(self.reference, sequence, self.match,
                                       self.mismatch, self.open, self.extend)
        (ref, aligned, score, mis, mat) = alignments[0]
        correlations = cor.correlate(ref, aligned)
        return {'reference': ref, 'sequence': aligned,
                'score': score, 'correlations': correlations}
