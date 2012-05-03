from subprocess import Popen, PIPE
import tempfile
import os
import re
import select

from rnastructure.secondary.connect import Parser as Connect
from rnastructure.secondary.dot_bracket import Parser as DotBracket


class FoldingTimeOutError(Exception):
    """This class indeciates that the folding operation took longer than the
    allowed time.
    """
    pass


class FoldingFailedError(Exception):
    """This class is used to indicate that the folding program exited with a
    non-zero status.
    """
    pass


class Folder(object):
    def __init__(self, directory=None, name='seq_file', length=40, time=120):
        self._filename = name
        self._base = directory
        self._time = time
        self._length = length
        self.names = []

    def generate_sequence_file(self, seq_file, sequence):
        """This method generates the file of sequences to be read by the
        folding programs.
        """
        seq_file.write(">sequence\n%s\n" % sequence)

    def generate_results(self, process, temp_dir, filename):
        """Generate the list or list like object which stores the results.
        """
        return ResultSet(temp_dir, filename)

    def validate_sequence(self, sequence):
        if not sequence:
            raise ValueError("Must give sequence to fold.")

        if len(sequence) > self._length:
            raise ValueError("Given sequence too long. Given: %s, Max: %s" %
                             (len(sequence), self._length))

    def fold(self, sequence, options=None):

        # Check the sequences are valid
        self.validate_sequence(sequence)

        temp_dir = self._base
        options = options or {}
        if not temp_dir:
            temp_dir = tempfile.mkdtemp()
        cur_dir = os.getcwd()

        # Write sequence file
        os.chdir(temp_dir)
        seq_file = open(self._filename, 'w')
        self.generate_sequence_file(seq_file, sequence)
        seq_file.close()

        # Run program
        args = [self.program]
        args.extend(self.process_arguments(self._filename, options))
        process = Popen(args, stdout=PIPE, stderr=PIPE)
        self.stdout = process.stdout
        self.stderr = process.stderr

        # Use select to wait until process is done.
        rlist, wlist, xlist = select.select([process.stderr], [],
                                            [process.stdout, process.stderr],
                                            self._time)
        os.chdir(cur_dir)

        if not rlist and not wlist and not xlist:
            process.kill()
            raise FoldingTimeOutError("Folding using %s timed out" %
                                      self.program)

        process.poll()
        code = process.returncode
        if code:
            raise FoldingFailedError("Fold program: %s failed. Status: %s " %
                                     (self.program, code))

        # Generate result to return.
        results = self.generate_results(process, temp_dir, self._filename)
        return results


class RNAalifold(Folder):
    """Use RNAalifold to fold some sequences. RNAalifold takes a sequence
    alignment and folds it to produce a secondary structure. For details on
    this program see: http://www.tbi.univie.ac.at/~ivo/RNA/ and
    http://www.tbi.univie.ac.at/~ronny/RNA/RNAalifold.html in particular.
    """
    program = 'RNAalifold'

    def generate_results(self, process, temp_dir, filename):
        """This will generate a list of size 1 because RNAalifold only
        generates a single structure. The result will be a Dot-Bracket parser
        with a sequence property set the consensus produced by RNAalifold and
        an energy property of the energy line given by RNAalifold.
        """

        self.raw = process.stdout.readlines()
        if len(self.raw) != 3 and len(self.raw) != 2:
            raise FoldingFailedError("No valid output")
        conensus = self.raw[0].rstrip()
        structure_line = self.raw[1].rstrip()
        parts = structure_line.split(' ', 1)
        parser = DotBracket(parts[0])
        parser.sequence = conensus
        parser.energy = parts[1]
        return [parser]

    def validate_sequence(self, sequences):
        if not sequences:
            raise ValueError("Must give sequence(s) to fold.")
        if isinstance(sequences, list):
            length = len(sequences[0])
            if length > self._length:
                raise ValueError("Given sequence too long. Given: %s, Max: %s"
                                 % (length, self._length))
            for sequence in sequences:
                if len(sequence) != length:
                    raise ValueError("All sequences must have the same length")
        else:
            super(RNAalifold, self).validate_sequence(sequences)

    def generate_sequence_file(self, seq_file, sequences):
        """This generates a stockholm like format that RNAalifold can read. It
        puts each sequence on a single line, which seems to be acceptable.
        """
        names = self.names
        if not names:
            names = ['sequence-%s' for index in xrange(len(sequences))]
        seq_file.write("# STOCKHOLM 1.0\n")
        for index, sequence in enumerate(sequences):
            name = names[index]
            seq_file.write("%s %s\n" % (name, sequence))

    def process_arguments(self, filename, options):
        args = []
        for key, value in options.items():
            prefix = '--'
            if len(key) == 1:
                prefix = '-'
            if value == True:
                args.append("%s%s" % (prefix, key))
            else:
                args.append("%s%s %s" % (prefix, key))
        args.append(filename)
        return args


class UNAfold(Folder):
    program = 'UNAFold.pl'

    def process_arguments(self, filename, options):
        args = []
        for key, value in options.items():
            prefix = '--'
            if len(key) == 1:
                prefix = '-'
            args.append('%s%s %s' % (prefix, key, value))
        args.append(filename)
        return args


class Mfold(Folder):
    program = 'mfold'

    def process_arguments(self, filename, options):
        args = []
        options['seq'] = filename
        for key, value in options.items():
            args.append("%s=%s" % (key.upper(), value))
        return args


class ResultSet(object):
    def __init__(self, base, name):
        self._name = name
        self._dir = base
        files = os.listdir(self._dir)
        pattern = '\A%s(_\d+)*\.ct\Z' % self._name
        valid = [name for name in files if re.match(pattern, name)]
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

    def __len__(self):
        return self._count


class Result(object):
    def __init__(self, name):
        self._name = name + ".%s"
        self.parser = Connect(self.connect_file())
        self.sequence = self.parser.sequence

    def connect_file(self):
        f = self.__file('ct')
        lines = f.readlines()
        f.close()
        return lines

    def indices(self, flanking=False):
        return self.parser.indices(flanking=flanking)

    def loops(self, flanking=False):
        return self.parser.loops(self.sequence, flanking=flanking)

    def __file(self, extension):
        ext_file = self._name % extension
        return open(ext_file, 'r')
