"""This is for storing utilities to wrap command line programs.
"""

import os
import sets
import select
import tempfile
from subprocess import Popen, PIPE


def is_number(arg):
    """Check that the given argument is a number.
    """
    return False


def is_true(arg):
    """Check that the given argument is true
    """
    return arg == True


class UnknownProgramOptionError(Exception):
    """This class represents given an unknown option to a program.
    """
    pass


class InvalidProgramOptionError(Exception):
    """This class represents a known option with an invalid value being given
    to the program.
    """
    pass


class InvalidInputError(Exception):
    """This class represents being given an invalid input.
    """


class ProgramTimeOutError(Exception):
    """This class indeciates that the program took longer than the allowed
    time.
    """
    pass


class ProgramFailedError(Exception):
    """This class is used to indicate that the program exited with a non-zero
    status.
    """
    pass


class Wrapper(object):
    """This is a class to ease wrapper of various command line programs.
    Classes inherting from this must define a program to run.
    """

    program = None

    options = {}

    def __init__(self, filename, directory=None, time=120):
        """Generate a new Wrapper.

        :filename: Filename to write to for input.
        :directory: Directory to write to. If not given then a tempdir is used.
        :time: Maximum time for the program.
        """
        self._filename = filename
        self._base = directory
        self._time = time
        self.stdout = None
        self.stderr = None
        self.known_options = set(self.options.keys())

    def _generate_input_file_(self, input_file, raw):
        """Generate the input file to process.
        """
        pass

    def _generate_results_(self, process, temp_dir, filename):
        """Generate a result object to use.
        """
        pass

    def validate(self, raw, options):
        """This validates the input and the options given.
        """
        if not self._validate_input_(raw):
            raise InvalidInputError("Input did not validate")
        if not self._validate_options_(options):
            raise InvalidProgramOptionError("Options did not validate")
        return True

    def _validate_input_(self, raw):
        """Validate the given raw input.
        """
        return True

    def _validate_options_(self, options):
        """Validate the given options.
        """
        extra_keys = []
        for key in options.keys():
            if key not in self.known_options:
                extra_keys.append(key)
        if extra_keys:
            key_str = "Options: %s are unknown." % ', '.join(extra_keys)
            raise UnknownProgramOptionError(key_str)

        for key, pattern in self.options.items():
            if key in options:
                value = options[key]
                if not pattern(value):
                    msg = "Key %s has invalid value %s" % (key, value)
                    raise InvalidProgramOptionError(msg)
        return True

    def _generate_arguments_(self, filename, options):
        """Process arguments to create the input argument array. This generates
        options in a very simple manner. If the option is only one character it
        generates a -opt followed by a value if necessary. If the option is long
        it generates --opt. Values are necessary if they are not True.
        """
        opts = []
        for key, value in options.iteritems():
            opt = '--%s'
            if len(key) == 1:
                opt = '-%s'
            opts.append(opt % value)
            if not is_true(value):
                opts.append(value)
        return opts

    def _program_failed_(self, process):
        return False

    def __call__(self, raw, options=None):
        """Run the program after generating the required input and return the
        produced results.
        """

        # Run program
        if not self.program:
            raise ValueError("Must define a program to run.")

        # Check the input and options are valid
        options = options or {}
        if not self.validate(raw, options):
            raise ValueError("Could not validate options and input")

        # Generate temp directory if needed.
        temp_dir = self._base
        if not temp_dir:
            temp_dir = tempfile.mkdtemp()
        cur_dir = os.getcwd()

        # Write input file
        os.chdir(temp_dir)
        with open(self._filename, 'w') as input_file:
            self._generate_input_file_(input_file, raw)

        args = [self.program]
        args.extend(self._generate_arguments_(self._filename, options))
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
            raise ProgramTimeOutError("Program %s timed out" % self.program)

        process.poll()
        code = process.returncode
        error_message = self._program_failed_(process)
        if code or error_message:
            msg = "Program %s failed. Status: %s. Message: %s." 
            raise ProgramFailedError(msg % (self.program, code, error_message))

        # Generate result to return.
        return self._generate_results_(process, temp_dir, self._filename)
