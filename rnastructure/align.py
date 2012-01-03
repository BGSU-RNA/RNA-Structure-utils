from suds import null
from suds.client import Client
from Bio import pairwise2 as nw


class Align(object):
    """Class to align RNA sequences. """

    def __init__(self, reference=None):
        # self._client = Client("http://www.rcsb.org/pdb/services/pdbws?wsdl")
        self.reference = reference
        self.match = 2
        self.mismatch = -1
        self.open = -0.5
        self.extend = -0.1

    def __lookup(self, sequence):
        if isinstance(sequence, dict):
            pdb = sequence['pdb']
            chain = sequence['chain']
            seq = client.service.getSequenceForStructureAndChain(pdb, chain)
            return seq
        return sequence

    def align_all(sequences):
        map(self.align, sequences)

    def best_alignment(self, sequences):
        alignments = self.align_all(sequences)
        return max(alignments, key=lambda a: a['score'])

    def align(self, sequence):
        seq = self.__lookup(sequence)
        alignments = nw.align.globalms(self.reference, seq, self.match,
                                       self.mismatch, self.open, self.extend)
        (ref, aligned, score, mis, mat) = alignments[0]
        correlations = self.correlate(ref, aligned)
        return {'reference': ref, 'sequence': aligned,
                'score': score, 'correlations': correlations}

    def correlate(self, reference, sequence):
        correlations = []
        for (i, char) in enumerate(reference):
            seq_char = sequence[i]
            if char != '-' and seq_char != '-':
                ref_pos = i - reference[:i].count('-')
                seq_pos = i - sequence[:i].count('-')
                correlations.append((ref_pos, seq_pos))
        return correlations
