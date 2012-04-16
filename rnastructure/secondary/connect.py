import re

import rnastructure.secondary.basic as basic


class Parser(basic.Parser):
    def __init__(self, lines):
        self.sequence = []
        pairs = self.__pairs(lines)
        self.sequence = ''.join(self.sequence)
        super(Parser, self).__init__(pairs)

    def __pairs(self, lines):
        pairs = []
        pattern = "\A(\d+)\s+([A-Za-z]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*\Z"
        for line in lines:
            if re.match(pattern, line):
                parts = line.split("\t")
                end = int(parts[4]) - 1
                if end < 0:
                    end = None
                # TODO: Assumes file is always sorted, is it?
                pairs.append(end)
                self.sequence.append(parts[1])
        return pairs
