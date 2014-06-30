from __future__ import with_statement

from os import path

from rnastructure.util.wrapper import is_true
from rnastructure.util.wrapper import is_number
from rnastructure.util.wrapper import Wrapper as Base
from rnastructure.secondary.basic import Parser as BaseParser
from rnastructure.secondary.connect import Parser as CtParser
from rnastructure.secondary.connect import Writer as CtWriter


class RemovePseudoknots(Base):
    """This is a class which wraps RemovePseudoknots. It is intended to remove
    the pseudoknots of an already parsed structure and return a newly parsed
    structure without the pseudoknots. See the .remove method for details.
    """
    program = "RemovePseudoknots"

    options = {
        'd': is_true,
        'D': is_true,
        'DNA': is_true,
        't': is_number,
        'T': is_number,
        'temperature': is_number
    }

    def __init__(self, directory=None, timeout=10):
        """Create a new Remove. This will remove pseudoknots from some RNA
        structure.

        :directory: Directory to work in. If none is given then a temp
        directory will be used.
        :timeout: Length of time to wait before the process times out.
        """
        self._output = 'output.ct'
        super(RemovePseudoknots, self).__init__('input.ct',
                                                directory=directory,
                                                time=timeout)

    def validate_input(self, raw):
        """Ensures that the input is a parsed secondary structure object.
        """
        return isinstance(raw, BaseParser)

    def generate_arguments(self, filename, options):
        """Creates the array of arguments to pass to RemovePseudoknots.
        """
        return [filename, self._output]

    def input_file(self, input_file, raw):
        """Takes the parsed secondary structure and generates the Connect file
        needed for RemovePseudoknots
        """
        formatter = CtWriter()
        formatter.write(input_file, raw)

    def results(self, process, temp_dir, filename):
        """Opens the result of RemovePseudoknots and creates a parsed
        secondary structure from it.
        """
        with open(path.join(temp_dir, self._output), 'r') as out:
            return CtParser(out)

    def _program_failed_(self, process):
        lines = process.stderr.readlines()
        if lines:
            del lines[0]  # First entry is a new line, so skip it.
            return ''.join(lines)
        return False

    def __call__(self, structure, options=None):
        """Remove pseudoknots from the given structure. This will return the
        resulting parsed structure. Note that RemovePseudoknots requires the
        DATADIR enviroment variable be set prior to running.

        :structure: The structure to remove pseudoknots from. This should be
        a parsed secondary structure like bpseq.Parser. It may be any one of
        those parser.
        :options: A hash of options to use.
        """
        return super(RemovePseudoknots, self).__call__(structure, options)
