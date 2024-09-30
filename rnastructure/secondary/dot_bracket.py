from collections import defaultdict
import string
import sys

import rnastructure.secondary.basic as basic


class Dialect(object):
    def __init__(self, unpaired, open_pair, close_pair, open_knot, close_knot):
        self.open_pair = open_pair
        self.close_pair = close_pair
        self.open_knot = open_knot
        self.close_knot = close_knot
        self.unpaired = unpaired

    def is_open_pair(self, char):
        return char in self.open_pair

    def is_close_pair(self, char):
        return char in self.close_pair

    def is_open_knot(self, char):
        return char in self.open_knot

    def is_close_knot(self, char):
        return char in self.close_knot

    def knot_type(self, char):
        if self.is_open_knot(char):
            return char
        elif self.is_close_knot(char):
            index = self.close_knot.index(char)
            open_index = len(self.close_knot) - 1 - index
            return self.open_knot[open_index]
        raise ValueError("Can't give knot_type of not a knot")

    def is_unpaired(self, char):
        return char in self.unpaired


class Parser(basic.Parser):
    """A class to parse 2D structures in dot-bracket format.

    The dot bracket format must be composed of <, (, [, {, ., :, }, ], ), >
    only. Paired bases are represented by (, <, >, and ). Unpaired bases are
    represented by ., and : while pseudonoted characters are represented by {,
    [, ], and }.
    """

    if sys.version_info[0] < 3:
        dialects = {
            'simple': Dialect('.:-', '(<', '>)', '{[', ']}'),
            'generic': Dialect('.:-', '(<', '>)', '{[' + string.uppercase,
                               string.lowercase[::-1] + ']}'),
            'rfam': Dialect('.;', '(<[{', '}]>)', string.uppercase,
                            string.lowercase[::-1]),
        }
    else:
        dialects = {
            'simple': Dialect('.:-', '(<', '>)', '{[', ']}'),
            'generic': Dialect('.:-', '(<', '>)', '{[' + string.ascii_uppercase,
                               string.ascii_lowercase[::-1] + ']}'),
            'rfam': Dialect('.;', '(<[{', '}]>)', string.ascii_uppercase,
                            string.ascii_lowercase[::-1]),
        }

    def __init__(self, structure, dialect='generic'):
        """Construct a new Parser object with the given structure.

        The structure should be a string in dot bracket notation with possible
        pseudoknots.

        """
        if isinstance(dialect, Dialect):
            self.__dialect = dialect
        elif dialect in self.dialects:
            self.__dialect = self.dialects[dialect]
        else:
            raise ValueError("Unknown dialect given")
        pairs = self.__pairs__(structure)
        super(Parser, self).__init__(pairs)

    def __getattr__(self, attr):
        if getattr(self.__dialect, attr):
            return getattr(self.__dialect, attr)
        return super(Parser, self).__getattr(attr)

    def __pairs__(self, structure):
        """Compute which bases are paired in the 2D structure.
        """
        helix_stack = []
        knot_stacks = defaultdict(list)
        pairs = [None] * len(structure)
        for index, char in enumerate(structure):
            if self.is_open_pair(char):
                helix_stack.append(index)
            elif self.is_close_pair(char):
                left = helix_stack.pop()
                pairs[left] = index
                pairs[index] = left
            elif self.is_unpaired(char):
                pass
            elif self.is_open_knot(char):
                knot_type = self.knot_type(char)
                knot_stacks[knot_type].append(index)
            elif self.is_close_knot(char):
                knot_type = self.knot_type(char)
                left = knot_stacks[knot_type].pop()
                pairs[left] = index
                pairs[index] = left
            else:
                raise ValueError("Unknown character: '%s'" % char)
        return pairs


class Writer(basic.Writer):
    def format(self, parser):
        dot_string = [None] * len(parser._pairs)
        for index, pair in enumerate(parser._pairs):
            if dot_string[index]:
                pass
            elif pair is None:
                dot_string[index] = '.'
            elif index < pair:
                dot_string[index] = '('
                dot_string[pair] = ')'
        return ''.join(dot_string)
