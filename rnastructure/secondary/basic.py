class EmptyStructureError(Exception):
    """This is a exception used with asked to parse something which has no
    pairs.
    """
    pass


class Writer(object):
    """Base class to format a parser structure as a string.
    """
    def write(self, open_file, parser):
        """Write the parser to a file.

        :open_file: The open file handle to write to.
        :parser: The parser to format.
        """
        return open_file.write(self.format(parser))

    def format(self, parser):
        """Create a string representation of the parser.

        :parser: The parser to format.
        """
        return str(parser._pairs)


class Parser(object):
    """This is the most generic parser for secondary structure. This builds
    with a list that gives the pairing information. This implements the actual
    algorithm for extracting loops. At the moment it can only extract hairpin
    and internal loops, in non-pseudoknotted structures but there is no reason
    why this cannot be extended to pseudoknotted structures and junctions.
    """

    def __init__(self, pairs, sequence=None):
        if not pairs:
            raise EmptyStructureError("Must specify pairs to find loops.")
        self.energy = ''
        self.sequence = sequence or [None] * len(pairs)
        self._pairs = pairs
        self._tree = Node((None, len(pairs)))
        self._loops = {}
        self.__as_tree()
        self.__find_indices(self._tree)

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

    def __find_indices(self, node):
        loop_type = node.loop_type()
        if loop_type:
            if loop_type not in self._loops:
                self._loops[loop_type] = []
            self._loops[loop_type].append(node.unpaired())
        for child in node.children:
            self.__find_indices(child)

    def loops(self, sequence=None, flanking=False):
        """Extract the loops for a given sequence. If no sequence is given then
        we try to use the sequence property of self, otherwise it is an error.

        :sequence: Sequence to extract loops from.
        :flanking: True if we wish to extract the flanking basepairs as well as
        the loop.
        """

        sequence = sequence or self.sequence
        if not sequence:
            raise ValueError("Must specify a sequence")

        if len(self) != len(sequence):
            msg = "Sequence has wrong size, given '%s' expected '%s'"
            raise ValueError(msg % (len(sequence), len(self)))

        def seq(parts, join_str='*'):
            if isinstance(parts[0], list):
                return join_str.join(map(lambda p: seq(p, ''), parts))
            return join_str.join(map(lambda p: sequence[p], parts))

        ranges = {}
        loops = self.indices(flanking=flanking)
        for name, total in loops.iteritems():
            ranges[name] = []
            for positions in total:
                char = '*'
                if name == 'hairpins':
                    char = ''
                loop_sequence = seq(positions, char)
                ranges[name].append(loop_sequence)
        return ranges

    def paired_base(self, index):
        """Get the base paired with the given one. None if no pair is made.
        """
        return self._pairs[index]

    def __flanking(self, part):
        """Get the flanking indices for the given part.
        """
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

    def __internal_flanking(self, loop):
        """Compute the flanking base pairs for the given internal loop.
        """
        left = self.__flanking(loop[0])
        right = []
        if len(loop) > 1:
            right = self.__flanking(loop[1])
        if not left:
            left = [self.paired_base(right[-1]), self.paired_base(right[0])]
        if not right:
            right = [self.paired_base(left[-1]), self.paired_base(left[0])]
        if left[0] > right[0]:
            return (right, left)
        return (left, right)

    def indices(self, flanking=False):
        """Get the indices of the loops in the parsed structure.

        :flanking: True if we wish to extract the positions of the flanking
        pairs as well.
        """

        if not flanking:
            return self._loops

        all_loops = {}
        for name, loops in self._loops.items():
            all_loops[name] = []
            for loop in loops:
                flank = None
                if name == "hairpins":
                    flank = self.__flanking(loop)
                elif name == 'internal':
                    flank = self.__internal_flanking(loop)
                else:
                    flank = tuple([self.__flanking(l) for l in loop])
                all_loops[name].append(flank)
        return all_loops

    def __len__(self):
        return len(self._pairs)


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
            if self.parent is None:
                return 'external'
            return 'internal'
        if len(self.children) > 1:
            return None
            # return 'junction'
        raise ValueError("Unknown type of loop")

    def left(self):
        if self.value[0] is None:
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
            child.print_tree(indent=indent + 1)

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
        return isinstance(other, Node) and self.value == other.value and \
            self.children == other.children
