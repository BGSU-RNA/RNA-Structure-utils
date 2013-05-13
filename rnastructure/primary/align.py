"""
This is a package for working with sequence alignments. This provides some
useful functions that we need for sequence alignments. We do extra work in this
package to make sure the alignments are unambiguous.
"""

from collections import namedtuple

from Bio.pairwise2 import align

from rnastructure.utils import correlate_aligned_indices as correlate


class AmbiguousAlignment(Exception):
    """This is raised when more than one possible alignment is detected and we
    cannot choose which is correct.
    """
    pass


class AlignmentFailed(Exception):
    """This is raised when we attempt an alignment but cannot perform one.
    """
    pass


Aligned = namedtuple('Aligned', ['reference', 'sequence', 'method', 'score'])


class Alignment(object):
    """A generic Alignment object. This object serves as a way to align a
    reference sequence and a target sequence. This will also generate the
    correlations between the names of positions in the reference and the
    target. If no names are given then we use the indexes.
    """

    def __init__(self, reference, sequence):
        """Create a new Alignment.

        :reference: The reference sequence.
        :sequence: The target sequence.
        """
        self._reference = reference
        self._sequence = sequence
        self.alignment = None
        self.correlations = {}

    def __call__(self):
        """Compute the alignment. This will actually execute the alignment
        between the reference and the sequence. This will add update the
        alignment property with the resulting alignment. In addition, this will
        update the correlations property with the correlations between the
        reference and sequence. If the reference and sequence have a names
        property they will be used to generate the correlations, otherwise just
        the indexes in the sequence will be used.
        """
        if self.alignment:
            return None
        self.alignment = self.__align__(self._reference, self._sequence)
        ref_names = getattr(self._reference, 'names', None)
        seq_names = getattr(self._sequence, 'names', None)
        self.correlations = correlate(self.alignment.reference,
                                      self.alignment.sequence,
                                      reference_names=ref_names,
                                      sequence_names=seq_names)


class NoGapsLocal(Alignment):
    """This creates a local alignment between the reference and the sequence
    where there are no gaps in the aligned region. We do allow mismatches,
    but never any gaps. When aligning this will first attempt an exact
    subsequence match. If we find more than one match this raises an
    AmbiguousAlignment exception. If we cannot find an exact match we then
    attempt to do a local alignment using a very high gap open and extend
    penalty. If we find more than one alignment with the optimal score we raise
    an AmbiguousAlignment exception. We also ensure that the alignment selected
    aligns at least the cutoff fraction of the target sequence. If the above
    methods fail we raise an AlignmentFailed exception.
    """

    def __init__(self, ref, seq, cutoff=0.9):
        """Create a new NoGapsLocal object.

        :ref: Reference sequence.
        :seq: Target sequence.
        :cutoff: The minimum fraction of the sequence which must be aligned.
        """

        self.cutoff = cutoff
        super(NoGapsLocal, self).__init__(ref, seq)

    def __make_aligned__(self, seq, start, stop):
        return '-' * start + seq[start:stop] + '-' * (len(seq) - stop)

    def __exact__(self, ref, seq):
        """Attempt an exact sequence match. This will only succeed if there is
        an exact substring in the reference sequence that is the target
        sequence. If there is more than one such substring an
        AmbiguousAlignment exception is raised.
        """

        start = ref.find(seq)
        if start != -1:
            if ref.rfind(seq) != start:
                raise AmbiguousAlignment("Multiple subsequence matches")
            overhang = len(ref) - len(seq)
            aligned_ref = self.__make_aligned__(ref, start, 0)
            aligned_seq = self.__make_aligned__(seq, start, overhang)
            return Aligned(aligned_ref, aligned_seq, 'exact', 1.0)
        return None

    def __gapless__(self, ref, seq):
        """Perform a gapless match. A gapless match is a local alignment with a
        very high gap open and extend penalty. If more than one local alignment
        with the maximal score is possible we raise an AmbiguousAlignment
        exception. If for some reason we have gaps despite the high penalty we
        return None.

        :ref: Reference sequence.
        :seq: Target sequence.
        """

        # -1000 is a very big penalty for opening or extending a gap, this
        # should prevent any gaps from being added. Though I suppose if the
        # sequence is very very long it may be possible.
        alignment = align.localxs(ref, seq, -1000, -1000)
        aligned_ref, aligned_seq, score, start, stop = alignment[0]

        maxes = filter(lambda x: x[2] == score, alignment)
        if len(maxes) > 1:
            raise AmbiguousAlignment("More than one optimal local alignment.")

        if stop - start != len(seq):
            return None

        count = 0.0
        for index, ref_char in enumerate(aligned_ref):
            if ref_char == aligned_seq[index]:
                count += 1.0

        return Aligned(alignment[0], alignment[1], 'gapless', count/len(seq))

    def __align__(self, ref, seq):
        """Get the alignment function for this object. In this case it is a
        function which first tries the exact method and then the gapless
        method. If both fail then a AlignmentFailed exceptions is raised.

        :ref: Reference sequence.
        :seq: Target sequence.
        """
        types = ['__exact__', '__gapless__']
        for name in types:
            func = getattr(self, name)
            result = func(ref, seq)
            if result and result.score > self.cutoff:
                return result
        raise AlignmentFailed("Could not generate an alignment.")

    def last(self):
        """Get the last aligned index in the reference sequence.
        """
        last = self.correlations[-1]
        if hasattr(self._reference, 'names'):
            return self._reference.names.index(last)
        return last


class GlobalFromLocal(Alignment):
    """This classes generates a global alignment out of a series of local
    alignments. We are given a series of short target sequences which are used
    to generate NoGapsLocal alignments to the reference sequence. We then use
    each of these to make a global alignment to the reference. Each small
    alignment cannot be overlapping. All target sequences must be completely
    aligned.
    """
    def __init__(self, first, second):
        """Create a new GlobalFromLocal alignment.
        The first should be the long sequence to align to. Second should be a
        list of short sequences to align.
        """
        super(GlobalFromLocal, self).__init__(first, second)
    pass
