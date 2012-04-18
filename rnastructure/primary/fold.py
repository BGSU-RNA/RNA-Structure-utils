from subprocess import Popen, PIPE
import tempfile
import os
import re
import select

from rnastructure.secondary.connect import Parser as Connect


class Folder(object):
    def __init__(self, directory=None, name='seq_file', time=120):
        self._filename = name
        self._base = directory
        self._time = time

    def fold(self, sequence, options=None):
        temp_dir = self._base
        options = options or {}
        if not temp_dir:
            temp_dir = tempfile.mkdtemp()
        cur_dir = os.getcwd()
        os.chdir(temp_dir)
        seq_file = open(self._filename, 'w')
        seq_file.write(">sequence\n%s\n" % sequence)
        seq_file.close()
        args = [self.program]
        args.extend(self.process_arguments(self._filename, options))
        process = Popen(args, stdout=PIPE, stderr=PIPE)
        rlist, wlist, xlist = select.select([process.stderr],
                                            [],
                                            [process.stdout, process.stderr],
                                            self._time)
        os.chdir(cur_dir)
        if not rlist and not wlist and not xlist:
            process.kill()
            raise ValueError("Failed running: %s" % self.program)
        return ResultSet(temp_dir, self._filename)


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
        self._sequence = None
        self._pairing = []

    def sequence(self):
        if not self._sequence:
            if not self._pairing:
                self.pairing()
            self._sequence = self._pairing.sequence
        return self._sequence

    def pairing(self):
        if not self._pairing:
            file_name = self._name % "ct"
            opened = open(file_name, 'r')
            self._pairing = Connect(opened)
            opened.close()

        return self._pairing
