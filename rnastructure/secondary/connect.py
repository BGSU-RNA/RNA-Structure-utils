"""This is a package for dealing with connect file formats of RNA secondary
structure. It contains code to write and read this format.
"""

import re

import rnastructure.secondary.basic as basic


class InvalidConnectLine(Exception):
    """This exceptions indicates that a line which is not a valid connect file
    line was found.
    """
    pass


class Writer(basic.Writer):
    """Format a parser as a connect file. I'm using the format of a connect
    file as described by:
        http://www.rnasoft.ca/strand/help.php
    and
        http://www.ibi.vu.nl/programs/k2nwww/static/data_formats.html

    A connect file has the following format. The first line is a header line
    of the form:
        $l Energy = $e
    where $l is the length of the parser and $e is the energy of the parser.
    The energy defaults to ''. The remaning lines are all of the form:
        $index $sequence $prev $next $pair $index
    where $index is the current 1 based index in the structure.

    $sequence is the sequence of the given index. This will default to '?',
    which may break some parsers. It is best to set a sequence property on the
    parser prior to writing to prevent this.

    $prev is the index of the previous position, $next is the index of the next
    position. If either index is less than 1 or greater than the length it is
    represented as 0.

    $pair is the index this position pairs with. Positions which do not pair
    are given as 0 here.
    """

    def format(self, parser):
        """Format the parser into a single string containing all the contents
        of a connect file.

        :parser: The parser to format.
        """
        header = '%s Energy = %s\n' % (len(parser), parser.energy)
        formatted = [header]
        for index in xrange(len(parser)):
            curr = index + 1
            after = curr + 1
            if curr >= len(parser):
                after = 0
            sequence = parser.sequence[index] or '?'
            pair = parser._pairs[index]
            if pair is None:
                pair = -1
            pair += 1
            data = (curr, sequence, index, after, pair, curr)
            line = '%s\t%s\t%s\t%s\t%s\t%s\n' % data
            formatted.append(line)
        return ''.join(formatted)


class Parser(basic.Parser):
    """Parse a connect file to get pairings. This will only take the first
    structure in the file.
    """

    def __init__(self, lines):
        self.sequence = []
        self.header = re.compile('\A(\d+)\s+(dG|Energy|ENERGY)')
        self.entry = \
            re.compile('\A(\d+)\s+([A-z?]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)')
        pairs = self.__pairs(lines)
        sequence = ''.join(self.sequence)
        super(Parser, self).__init__(pairs, sequence=sequence)

    def __pairs(self, lines):
        pairs = []
        for index, line in enumerate(lines):
            line = line.strip()
            if index == 0 and self.header.match(line):
                pass
            elif index > 0 and self.header.match(line):
                break
            elif self.entry.match(line):
                parts = line.split()
                end = int(parts[4]) - 1
                if end < 0:
                    end = None
                # TODO: Assumes file is always sorted, is it?
                pairs.append(end)
                self.sequence.append(parts[1])
            else:
                raise InvalidConnectLine("Invalid line: %s" % line)
        return pairs
