import re

import rnastructure.secondary.basic as basic


class Parser(basic.Parser):
    """Parse a connect file to get pairings. This will only take the first
    structure in the file.
    """

    def __init__(self, lines):
        self.sequence = []
        pairs = self.__pairs(lines)
        if not pairs:
            raise ValueError, "Could not parse given data"
        self.sequence = ''.join(self.sequence)
        super(Parser, self).__init__(pairs)

    def __pairs(self, lines):
        pairs = []
        pattern = r'\A(\d+)\s+([A-Za-z]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*'
        end_pattern = r'\A\d+\s+dG\s*=\s*'
        for index,line in enumerate(lines):
            if re.match(end_pattern, line) and index > 0:
                break
            if re.match(pattern, line):
                parts = line.split("\t")
                end = int(parts[4]) - 1
                if end < 0:
                    end = None
                # TODO: Assumes file is always sorted, is it?
                pairs.append(end)
                self.sequence.append(parts[1])
        return pairs
