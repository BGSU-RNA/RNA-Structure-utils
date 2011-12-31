from itertools import takewhile


class Parser(object):
    """
    A Simple class to parse 2D structures.
    This cannot handle pseudoknots.
    """

    def __init__(self, structure):
        self.loops = {'hairpins': [], 'internal': []}
        self._len = len(structure)

        print(structure)
        pairs = self.__map_relations(structure)
        print(pairs)
        node = self.__convert(pairs, Node())
        print(node)
        self.__find_loops(node)
        print(self.loops)

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

        print('pairs', pairs)
        start = 0
        if left:
            start = left[-1][0] + 1
        if start < len(pairs):
            end = start + pairs[start][0]
            node.add(self.__convert(pairs[(start + 1):end], Node()))
            while start <= end:
            # for i in range(pairs[start][0]):
                start = start + 1

        print('start', start)
        print('rest', pairs[start:])

        self.__convert(pairs[start:], node)
        # [node.add(r[-1]) for r in pairs[start:] if not r[0]]

        return node

    def __find_loops(self, node):
        loop = node.get_loop()
        if len(loop) == 1:
            print('hairpin', loop)
            self.loops['hairpins'].append(loop[0])
        else:
            if len(loop) == 2:
                print('internal', loop)
                self.loops['internal'].append(loop)
            [self.__find_loops(n) for n in node if isinstance(n, Node)]

    def parse(self, sequence):
        if len(self) != len(sequence):
            msg = "Bad sequence length of dotbracket, given {0} expected {1}"
            raise ValueError(msg.format(len(sequence), len(self)))

        ranges = {}
        for name, total in self.loops.iteritems():
            ranges[name] = []
            for positions in total:
                seq = ''.join(map(lambda v: sequence[v], positions))
                ranges[name].append(seq)
        return ranges

    def __len__(self):
        return self._len


class Node(object):
    def __init__(self):
        self._components = []

    def add(self, value):
        self._components.append(value)

    def get_loop(self):
        def take_loop(it):
            return takewhile(lambda c: not isinstance(c, Node), it)
        left = list(take_loop(self))

        if left == self._components:
            return [tuple(self._components)]
        else:
            right = list(take_loop(reversed(self._components)))
            if right:
                right.reverse()
                return tuple([left, right])

        return ()

    def __iter__(self):
        return iter(self._components)

    def __len__(self):
        return len(self._components)

    def __getitem__(self, index):
        return self._components[index]
