from itertools import takewhile


class Parser(object):
    """
    A Simple class to parse 2D structures.
    This cannot handle pseudoknots.
    """

    def __init__(self, structure):
        self.loops = {'hairpins': [], 'internal': [], 'junction': [],
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
                raise ValueError("Unknown character: '{0}'".format(char))
        return pairs

    def __convert(self, pairs, node):
        if not pairs:
            return node
        left = takewhile(lambda (_, p): not p[0], enumerate(pairs))
        left = list(left)
        [node.append((l[1], l[-1])) for (_, l) in left]

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
            self.loops['hairpins'].append(loop[0])
        elif loop:
            if len(loop) == 2:
                self.loops['internal'].append(loop)
            else:
                self.loops['junction'].append(loop)

        knot = node.pseudoknot()
        if knot:
            self.loops['pseudoknot'].append(knot)
        [self.__find_loops(n) for n in node if isinstance(n, Node)]

    def parse(self, sequence):
        if len(self) != len(sequence):
            msg = "Bad sequence length of dotbracket, given {0} expected {1}"
            raise ValueError(msg.format(len(sequence), len(self)))

        def seq(parts, join_str='*'):
            if isinstance(parts[0], list):
                return join_str.join(map(lambda p: seq(p, ''), parts))
            return join_str.join(map(lambda p: sequence[p], parts))

        ranges = {}
        for name, total in self.loops.iteritems():
            ranges[name] = []
            for positions in total:
                ch = '*'
                if name == 'hairpins':
                    ch = ''
                loop_sequence = seq(positions, ch)
                ranges[name].append(loop_sequence)
        return ranges

    def __len__(self):
        return self._len


class Node(list):
    def __find(self, fn):
        result = [[]]
        for (i, entry) in enumerate(self):
            if isinstance(entry, Node):
                if result[-1]:
                    result.append([])
            else:
                if fn(entry):
                    result[-1].append(entry[1])
        if result[0] == []:
            return tuple()
        return tuple(result)

    def pseudoknot(self):
        return self.__find(lambda (c, i): c == '{' or c == '}')

    def loop(self):
        return self.__find(lambda (c, i): c == '.')
