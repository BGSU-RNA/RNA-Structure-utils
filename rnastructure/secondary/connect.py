import re

import rnastructure.secondary.basic as basic


class InvalidConnectLine(Exception):
    """This exceptions indicates that a line which is not a valid connect file
    line was found.
    """


class Parser(basic.Parser):
    """Parse a connect file to get pairings. This will only take the first
    structure in the file.
    """

    def __init__(self, lines):
        self.sequence = []
        self.header = re.compile('\A\d+\s+dG\s*=\s*')
        self.entry = \
          re.compile('\A(\d+)\s+([A-z]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*')
        pairs = self.__pairs(lines)
        self.sequence = ''.join(self.sequence)
        super(Parser, self).__init__(pairs)

    def __pairs(self, lines):
        pairs = []
        for index, line in enumerate(lines):
            if index == 0 and self.header.match(line):
                pass
            elif index > 0 and self.header.match(line):
                break
            elif self.entry.match(line):
                parts = line.split("\t")
                end = int(parts[4]) - 1
                if end < 0:
                    end = None
                # TODO: Assumes file is always sorted, is it?
                pairs.append(end)
                self.sequence.append(parts[1])
            else:
                raise InvalidConnectLine("Invalid line: %s" % line)
        return pairs
