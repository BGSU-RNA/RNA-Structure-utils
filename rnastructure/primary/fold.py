"""This module wraps up several common and useful folding programs. In general
these programs take one or more sequences and fold them into one or more
secondary structures. Each program is represented as a single class where
objects from each class are callable, take the appropriate inputs and return a
parsed secondary structure object from rnastructure.secondary.

Classes are named according to what program they wrap. We always capitalize the
first character of each name but otherwise all names are as the program names.
"""

from __future__ import with_statement


import os
import re
from collections import MutableSequence

from rnastructure.util.wrapper import Wrapper
from rnastructure.util.wrapper import InvalidInputError
from rnastructure.secondary.connect import Parser as Connect
from rnastructure.secondary.dot_bracket import Parser as DotBracket
from rnastructure.secondary.rnaplot import PostScriptParser as RNAPlot


class FoldingFailedError(Exception):
    pass


class Folder(Wrapper):
    """Folder wraps up several of the common patterns for running the folding
    programs. It is uses the Wrapper classes to give a uniform interface.
    """

    def __init__(self, length=500, directory=None, name='seq_file', time=120):
        """Create a new Folder.

        :length: Maximum length of the sequence to allow.
        :directory: Directory to work in, a temp dir if None is given.
        :name: Name of the input file to use.
        :time: Maximum amount of time to give the folding program.
        """
        self._length = length
        self.sequence_names = []
        super(Folder, self).__init__(name, directory=directory, time=time)

    def input_file(self, seq_file, sequence):
        """Generate a simple fasta like format for the sequences. This differs
        from a strict fasta format in that sequence lines are not wrapped at 80
        characters and there is no way to generate a header.

        :seq_file: File object to write to.
        :sequence: Sequence to write.
        """
        seq_file.write(">sequence\n%s\n" % sequence)

    def validate_input(self, sequence):
        """Check that the the sequence is a string, that it has the correct
        length and that it is only composed of A, C, G, U, and '-' ignoring
        case.

        :sequence: Sequence to check.
        """

        if len(sequence) > self._length:
            msg = "Given sequence too long. Given: %s, Max: %s"
            raise InvalidInputError(msg % (len(sequence), self._length))

        if not re.match("[acguACGU-]+", sequence):
            raise InvalidInputError("Sequence may only be A, C, G, U, -")

        return True

    def results(self, process, temp_dir, filename):
        """Generate a ResultSet object for the results.

        :process: The process object.
        :temp_dir: Directory all work was done in.
        :filename: Input filename.
        """
        return ResultSet(temp_dir, filename)


class RNAalifold(Folder):
    """Use RNAalifold to fold some sequences. RNAalifold takes a sequence
    alignment and folds it to produce a single secondary structure. For
    details on this program see: http://www.tbi.univie.ac.at/~ivo/RNA/ and
    http://www.tbi.univie.ac.at/~ronny/RNA/RNAalifold.html in particular.
    """

    program = 'RNAalifold'

    def results(self, process, temp_dir, filename):
        """This will generate a list of size 1 because RNAalifold only
        generates a single structure. The result will be a Dot-Bracket parser
        with a sequence property set the consensus produced by RNAalifold and
        an energy property of the energy line given by RNAalifold.

        :process: The process object.
        :temp_dir: The directory all work was done in.
        :filename: Input filename.
        """
        self.raw = process.stdout.readlines()
        if len(self.raw) != 3 and len(self.raw) != 2:
            raise FoldingFailedError("No valid output")
        consensus = self.raw[0].rstrip()
        structure_line = self.raw[1].rstrip()
        parts = structure_line.split(' ', 1)
        parser = DotBracket(parts[0])
        parser.sequence = consensus
        parser.energy = parts[1]

        ps_parser = self._load_locations(temp_dir, filename)
        parser.locations = ps_parser.locations
        parser.box = ps_parser.box

        return [parser]

    def _load_locations(self, temp_dir, filename):
        ps_file = os.path.join(temp_dir, 'alirna.ps')
        with open(ps_file, 'r') as raw:
            return RNAPlot(raw)

    def validate_input(self, sequences):
        """Check that the input is a list of valid sequences of the same
        length.

        :sequences: Input sequences to check.
        """
        if not isinstance(sequences, list) or len(sequences) < 1:
            raise InvalidInputError("Must give a list of sequences to fold.")
        length = len(sequences[0])
        for sequence in sequences:
            if len(sequence) != length:
                raise InvalidInputError("All sequences have the same length")
            super(RNAalifold, self).validate_input(sequence)
        return True

    def input_file(self, seq_file, sequences):
        """This generates a stockholm like format that RNAalifold can read. It
        puts each sequence on a single line, which seems to be acceptable. This
        only generates the required STOCKHOLM 1.0 header and then sequence
        lines. Each sequence line is of the form: 'sequence-$n $sequence'.
        Where $n is the index and $sequence is the sequence at that index.

        :seq_file: File to write to.
        :sequences: The input sequences.
        """
        names = self.sequence_names
        if not names:
            names = ['sequence-%s' for index in xrange(len(sequences))]
        seq_file.write("# STOCKHOLM 1.0\n")
        for index, sequence in enumerate(sequences):
            name = names[index]
            seq_file.write("%s %s\n" % (name, sequence))

    def __call__(self, sequences, options=None):
        """Fold the given sequences with RNAalifold. All sequences must have
        the same length and should be only A, C, G, U or -. There must be at
        least 2 sequences to fold. For details of what can be raised see
        rnastructure.util.wrapper.Wrapper.__call__.

        Options are currently ignored.

        :sequences: List of sequences to fold.
        :options: Options for RNAalifold.
        :returns: A parsed secondary structure object.
        """
        return super(RNAalifold, self).__call__(sequences)


class UNAFold(Folder):
    """This class wraps up UNAFold for use. UNAFold is the sucessor to mfold.
    UNAFold takes a single sequence and folds it to produce a secondary
    structure. For details see: http://mfold.rna.albany.edu/.
    """

    program = 'UNAFold.pl'

    def __call__(self, sequence, options=None):
        """Fold the given sequence with UNAFold.

        Currently options are ignored.

        :sequence: Sequence to fold.
        :options: The program options.
        :returns: A parsed secondary structure object.
        """
        return super(UNAFold, self).__call__(sequence)


class ResultSet(MutableSequence):
    def __init__(self, base, name):
        self._name = name
        self._dir = base
        files = os.listdir(self._dir)
        pattern = '%s(_\d+)*\.ct' % self._name
        valid = [n for n in files if re.match(pattern, n)]
        if not valid:
            raise FoldingFailedError("Could not find any generated foldings")
        self._count = len(valid)
        self._pairings = [None] * self._count

    def __getitem__(self, index):
        if index >= self._count:
            raise IndexError("Index out of bounds")
        if not self._pairings[index]:
            file_name = self._name
            if index > 0:
                file_name += "_%s" % index
            full_name = os.path.join(self._dir, file_name)
            self._pairings[index] = Result(full_name)
        return self._pairings[index]

    def __setitem__(self, index, value):
        self._pairings = value

    def __delitem__(self, index):
        del self._pairings[index]

    def insert(self, index, value):
        self._pairings.insert(index, value)

    def __len__(self):
        return self._count


class Result(object):
    def __init__(self, name):
        self._name = name + ".%s"
        self.parser = Connect(self.connect_file())
        self.sequence = self.parser.sequence

    def connect_file(self):
        with self.__file__('ct') as f:
            return f.readlines()

    def indices(self, flanking=False):
        return self.parser.indices(flanking=flanking)

    def loops(self, flanking=False):
        return self.parser.loops(self.sequence, flanking=flanking)

    def __file__(self, extension):
        ext_file = self._name % extension
        return open(ext_file, 'r')
