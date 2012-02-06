import re

import rnastructure.secondary.basic as basic


class Parser(basic.Parser):
    def __init__(self, lines):
        self.sequence = []
        pairs = self.__pairs(lines)
        self.sequence = ''.join(self.sequence)
        super(Parser, self).__init__(pairs)

    def __pairs(self, lines):
        pattern = "\A(\d+)\s+([a-zA-Z]+)\s+(\d+)\Z"
        pairs = []
        for line in lines:
            if re.match(line, pattern):
                parts = line.split("\t")
                self.sequence.append(parts[1])
                end = int(parts[2]) - 1
                if end < 0:
                    end = None
                pairs.append(end)
        return pairs
