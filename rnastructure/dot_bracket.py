from itertools import takewhile


class Parser(object):
    """
      A Simple class to parse 2D structures.
    """

    def __init__(self, structure):
        """
          Construct a new Parser object with the given structure. The structure
          should be a string in dot bracket notation with possible pseudoknots.
        """
        self._loops = {'hairpins': [], 'internal': [], 'junction': [],
                      'pseudoknot': []}
        self._len = len(structure)

        pairs = self.__map_relations(structure)
        node = self.__convert(pairs, Node())
        self.__find_loops(node)

    def __map_relations(self, structure):
        helix_stack = []
        knot_stack = []
        pairs = [False] * len(structure)
        for index, char in enumerate(structure):
            if char == '(':
                helix_stack.append(index)
            elif char == ')':
                left = helix_stack.pop()
                pairs[left] = (index - left, '(', left)
                pairs[index] = (left - index, ')', index)
            elif char == '.':
                pairs[index] = (False, '.', index)
            elif char == '{':
                knot_stack.append(index)
            elif char == '}':
                left = knot_stack.pop()
                pairs[left] = (False, '{', left)
                pairs[index] = (False, '}', index)
            else:
                raise ValueError("Unknown character: '%s'" % char)
        return pairs

    def __convert(self, pairs, node):
        if not pairs:
            return node
        left = takewhile(lambda (_, p): not p[0], enumerate(pairs))
        left = list(left)
        for (_, entry) in left:
            node.append((entry[1], entry[-1]))

        start = 0
        if left:
            start = left[-1][0] + 1
        if start < len(pairs):
            end = start + pairs[start][0]
            node.append(self.__convert(pairs[(start + 1):end], Node()))
            while start <= end:
                start = start + 1

        self.__convert(pairs[start:], node)

        return node

    def __find_loops(self, node):
        loop = node.loop()
        if loop and len(loop) == 1:
            self._loops['hairpins'].append(loop[0])
        elif loop:
            if len(loop) == 2:
                self._loops['internal'].append(loop)
            else:
                self._loops['junction'].append(loop)

        knot = node.pseudoknot()
        if knot:
            self._loops['pseudoknot'].append(knot)

        for entry in node:
            if isinstance(entry, Node):
                self.__find_loops(entry)

    def parse(self, sequence):
        if len(self) != len(sequence):
            msg = "Bad sequence length of dotbracket, given '%s' expected '%s'"
            raise ValueError(msg % (len(sequence), len(self)))

        def seq(parts, join_str='*'):
            if isinstance(parts[0], list):
                return join_str.join(map(lambda p: seq(p, ''), parts))
            return join_str.join(map(lambda p: sequence[p], parts))

        ranges = {}
        for name, total in self._loops.iteritems():
            ranges[name] = []
            for positions in total:
                char = '*'
                if name == 'hairpins':
                    char = ''
                loop_sequence = seq(positions, char)
                ranges[name].append(loop_sequence)
        return ranges

    def loops(self, flanking=False):
        def add_flank(part):
            if not part:
                return part
            flank = list(part)
            left = part[0]
            right = part[-1]
            if left > 0:
                flank.insert(0, left - 1)
            if right + 1 < len(self):
                flank.append(right + 1)
            return flank

        if not flanking:
            return self._loops

        all_loops = {}
        for name, loops in self._loops.items():
            all_loops[name] = []
            for loop in loops:
                if name == "hairpins":
                    all_loops[name].append(add_flank(loop))
                else:
                    flank = [add_flank(l) for l in loop]
                    all_loops[name].append(tuple(flank))
        return all_loops

    def __len__(self):
        return self._len


class Node(list):
    def __find(self, func):
        result = [[]]
        for (_, entry) in enumerate(self):
            if isinstance(entry, Node):
                if result[-1]:
                    result.append([])
            else:
                if func(entry):
                    result[-1].append(entry[1])
        if not result[0]:
            return tuple()
        return tuple(result)

    def pseudoknot(self):
        """
          Find the pseudoknot, if any, in this node. Pseudoknots entries in the
          node with either '{' or '}' characters. If a pseudoknot is found it
          will be returned as tuple of the form (left_half, right_half). If no
          pseudoknot is found then an empty tuple is returned.
        """
        return self.__find(lambda (c, i): c == '{' or c == '}')

    def loop(self):
        """
          Find the loop in this node if any. Nodes are entries in the node with
          '.' character. If a loop is found it will be returned as a tuple of
          the form (first, second, ...). If no loop is found an empty tuple is
          returned.
        """
        return self.__find(lambda (c, i): c == '.')
