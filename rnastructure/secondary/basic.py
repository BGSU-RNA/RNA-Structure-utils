

class Format(object):
    def __init__(self, parser):
        self._pairs = parser._pairs

    def __str__(self):
        return self.format()


class Parser(object):
    def __init__(self, pairs):
        if not pairs:
            raise ValueError("Must specify pairs to find loops.")
        self._pairs = pairs
        self._len = len(pairs)
        self._tree = Node((None, self._len))
        self._loops = {}
        self.__as_tree()
        self.__find_loops(self._tree)

    def __as_tree(self):
        stack = []
        for i, pair in enumerate(self._pairs):
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

    def __find_loops(self, node):
        loop_type = node.loop_type()
        if loop_type:
            if loop_type not in self._loops:
                self._loops[loop_type] = []
            self._loops[loop_type].append(node.unpaired())
        for child in node.children:
            self.__find_loops(child)

    def parse(self, sequence):
        if len(self) != len(sequence):
            msg = "Sequence has wrong size, given '%s' expected '%s'"
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
            return self.children[-1]
        return Node((None, None))

    def loop_type(self):
        if not self.unpaired():
            return None
        if not self.children:
            return 'hairpin'
        if len(self.children) == 1:
            return 'internal'
        if len(self.children) > 1:
            return 'junction'
        raise ValueError("Unknown type of loop")

    def left(self):
        if self.value[0] == None:
            return 0
        return self.value[0] + 1

    def spans(self, flanking=False):
        start = self.left()
        end = self.value[1]
        if flanking:
            start -= 1
            end += 1
        return range(start, end)

    def unpaired(self):
        unpaired = []
        left = self.left()
        for child in self.children:
            right = child.value[0]
            looped = range(left, right)
            if looped:
                unpaired.append(looped)
            left = child.value[1] + 1
        last = range(left, self.value[1])
        if last:
          unpaired.append(last)
        return tuple(unpaired)

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

    def print_tree(self, indent=0):
        print(" " * indent + "Node: " + str(self.value))
        for child in self.children:
          child.print_tree(indent=indent+1)

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
