from itertools import takewhile


class Parser(object):
    """
    A Simple class to parse 2D structures.
    This cannot handle pseudoknots.
    """

    def __init__(self, structure):
        self.loops = {'hairpins': [], 'internal': [], 'junction': []}
        self._len = len(structure)

        pairs = self.__map_relations(structure)
        node = self.__convert(pairs, Node())
        self.__find_loops(node)

    def __map_relations(self, structure):
        stack = []
        pairs = [False] * len(structure)
        for index, char in enumerate(structure):
            if char == '(':
                stack.append(index)
            elif char == ')':
                left = stack.pop()
                pairs[left] = (index - left, '(', left)
                pairs[index] = (left - index, ')', index)
            elif char == '.':
                pairs[index] = (False, '.', index)
            else:
                raise ValueError("Unknown character: '{0}'".format(char))
        return pairs

    def __convert(self, pairs, node):
        if not pairs:
            return node
        left = takewhile(lambda (_, p): not p[0], enumerate(pairs))
        left = list(left)
        [node.add(l[-1]) for (_, l) in left]

        start = 0
        if left:
            start = left[-1][0] + 1
        if start < len(pairs):
            end = start + pairs[start][0]
            node.add(self.__convert(pairs[(start + 1):end], Node()))
            while start <= end:
                start = start + 1

        self.__convert(pairs[start:], node)
        # [node.add(r[-1]) for r in pairs[start:] if not r[0]]

        return node

    def __find_loops(self, node):
        loop = node.get_loop()
        if loop and len(loop) == 1:
            self.loops['hairpins'].append(loop[0])
        elif loop:
            if len(loop) == 2:
                self.loops['internal'].append(loop)
            else:
                self.loops['junction'].append(loop)
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


class Node(object):
    def __init__(self):
        self._components = []

    def add(self, value):
        self._components.append(value)

    def get_loop(self):
        loop = [[]]
        for (i, entry) in enumerate(self):
            if isinstance(entry, Node):
                if loop[-1]:
                    loop.append([])
            else:
                loop[-1].append(entry)

        if loop[0] == []:
            return tuple()
        return tuple(loop)

    def __iter__(self):
        return iter(self._components)

    def __len__(self):
        return len(self._components)

    def __getitem__(self, index):
        return self._components[index]
