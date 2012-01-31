

class Basic(object):
    def __init__(self, pairs):
        self._pairs = pairs
        self._len = len(pairs)
        self._tree = Node((None, self._len))
        self._loops = {'hairpin': [],
                       'internal': [],
                       'junction': [],
                      'pseudoknot': []}
        self.__as_tree()
        self.__find_loops()

    def __as_tree(self):
        """Create the parse tree for this 2D structure.
        """
        stack = []
        for i in range(self._len):
            pair = self._pairs[i]
            if i < pair:
                stack.append((i, pair))
            elif None < pair < i:
                end = i
                while stack and stack[-1][0] > pair:
                    end = max(end, stack.pop()[1])
                stack[-1] = (stack[-1][0], max(end, stack[-1][1]))

            if stack and i == stack[-1][1]:
                pair = stack.pop()
                node = Node(pair)
                self._tree.add_to_tree(node)

    def __find_loops(self):
        pass

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


class Node(object):
    def __init__(self, value, parent=None):
        self.value = value
        self.parent = parent
        self.children = []

    def largest(self):
        if self.children:
            return self.children[0]
        return Node((None, None))

    def add_to_tree(self, child):
        biggest = self.largest()
        while biggest > child:
            child.add_child(biggest)
            self.children.remove(biggest)
            biggest = self.largest()
        self.add_child(child)

    def add_child(self, child):
        child.parent = self
        self.children.append(child)

    def __ne__(self, other):
        return self.value != other.value

    def __ge__(self, other):
        return self.value >= other.value

    def __gt__(self, other):
        return self.value > other.value

    def __le__(self, other):
        return self.value <= other.value

    def __lt__(self, other):
        return self.value < other.value

    def __eq__(self, other):
        return self.value == other.value and self.children == other.children
