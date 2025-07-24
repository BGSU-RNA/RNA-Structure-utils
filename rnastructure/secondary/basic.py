"""
Code to parse secondary structures in dot-bracket format, and extract sequences
"""

class EmptyStructureError(Exception):
    """This is an exception used with asked to parse something which has no
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
    """
    This is the most generic parser for secondary structure. This builds
    with a list that gives the pairing information. This implements the actual
    algorithm for extracting loops.
    It was originally written to extract hairpin and internal loops,
    in non-pseudoknotted structures.
    It was extended to J3 and J4 in 2025
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
        """
        Work through the list of pairing partner in self._pairs,
        identifying when a new branch is formed and adding that to a tree.

        self._pairs tells what each position is paired with, like:
        [18, None, None, None, 7, None, None, 4, None, None, None, 14, None, None, 11, None, None, None, 0]
        it is computed in dot_bracket before this is called
        """

        # print('\nbasic.py: %s' % self._pairs)

        stack = []
        for i, pair in enumerate(self._pairs):
            if pair is None:
                continue
            elif i < pair:
                # pair is opening
                stack.append((i, pair))
            elif pair < i:
                # pair is closing
                end = i
                while stack and stack[-1][0] > pair:
                    end = max(end, stack.pop()[1])
                stack[-1] = (stack[-1][0], max(end, stack[-1][1]))

            if stack and i == stack[-1][1]:
                pair = stack.pop()
                node = Node(pair)

                # print('__as_tree is popping pair %s,%s off the stack' % pair)
                # print('Tree from this node')
                # node.print_tree()

                self._tree.add_to_tree(node)

        # print('Here is the whole tree, which lists every basepair')
        # print('Basepairs are indented according to how deep they are in the secondary structure')
        # self._tree.print_tree()

    def __find_indices(self, node):
        """
        Starting at the given node, recursively go through children
        and identify the ranges of nucleotides in each loop.
        """
        loop_type = node.loop_type()
        if loop_type:
            if loop_type not in self._loops:
                self._loops[loop_type] = []
            # previously:
            # node.unpaired returns only the strands with non-empty interior
            # it will miss entire strands for J3, J4 in some cases!
            # new:
            # self._loops[loop_type] is a list of tuples of start and end positions
            # of the unpaired/interior bases in each strand
            self._loops[loop_type].append(node.unpaired())
        for child in node.children:
            self.__find_indices(child)

    def loops(self, sequence=None, flanking=False):
        """
        Extract the loops for a given sequence. If no sequence is given then
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

        type_to_loop_sequences = {}
        type_to_loops = self.indices(flanking=flanking)
        for name, loops in type_to_loops.items():
            # name is loop type, loops is a list of loops
            # print('basic.py: name: %s loops: %s' % (name,str(loops)))

            type_to_loop_sequences[name] = []
            for loop in loops:
                # print('basic.py: loop: %s' % str(loop))
                if name == 'hairpin':
                    char = ''
                else:
                    char = '*'
                loop_sequence = seq(loop, char)
                type_to_loop_sequences[name].append(loop_sequence)
        return type_to_loop_sequences

    def paired_base(self, index):
        """Get the base paired with the given one. None if no pair is made.
        """
        return self._pairs[index]

    def __flanking(self, part):
        """
        Get the flanking indices for the given part.
        Converts a range or list of unpaired bases on a single strand to be a little
        longer and include the paired bases as well.
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

    def __new_flanking(self, left, right):
        """
        Get the flanking indices for the given strand.
        left and right are the first and one past the last unpaired bases on the strand
        longer and include the paired bases as well.
        It can happen that left == right when there are no unpaired bases.
        In that case, use the value of right to get the correct indices.
        """

        # start with the unpaired nucleotides
        if left == right:
            flank = []
        else:
            flank = list(range(left,right))

        if left > 0:
            flank.insert(0, left - 1)
        if right < len(self):
            flank.append(right)

        return flank

    def __internal_flanking(self, loop):
        """
        Compute the flanking base pairs for the given internal loop.
        It may be that one of the strands is not present because it has length zero.
        This code works around that possibility for IL
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
        """
        Loops are already extracted and stored in self._loops as upper
        and lower indices of the unpaired nucleotides of each strand.
        This method fills in the indices of each strand.
        Go over each loop type, then over each loop, then over each strand
        When flanking is True, flanking pairs are added.

        :flanking: True if we wish to extract the positions of the flanking
        pairs as well.
        """

        # new data structure is more complicated, can't just return self._loops
        # if not flanking:
            # return self._loops

        all_loops = {}
        for name, loops in self._loops.items():
            # name is a loop type like HL, IL, J3, J4
            # loops is a list of all of the loops of that type

            # by this point, loop types and indices have already been determined
            # print("basic.py: name: %s unpaired ranges: %s" % (name, loops))

            all_loops[name] = []
            for loop in loops:
                # print('basic.py: loop_type: %s' % name)

                # old code needed special treatment for different cases
                # flank = None
                # if name == "hairpins":
                #     flank = self.__flanking(range(a,b))
                # elif name == 'internal':
                #     # special code for internal loops deals with empty ranges
                #     flank = self.__internal_flanking(range(a,b))
                # else:
                #     flank = tuple([self.__flanking(l) for l in loop])

                # new code has the information needed for each strand
                strands = []
                for a,b in loop:
                    if flanking:
                        # extend strand to flanking pairs
                        strands.append(self.__new_flanking(a,b))
                    else:
                        # include only the unpaired nucleotides; range may be empty
                        strands.append(range(a,b))

                # print('basic.py: loop:      %s' % str(loop))
                # print('basic.py: indices:   %s' % str(flank))

                if len(strands) > 0:
                    strands = tuple(strands)
                else:
                    strands = None

                all_loops[name].append(strands)
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
            # for example, successive Watson-Crick basepairs
            return None
        if self.parent is None:
            # outside the first Watson-Crick basepair, or
            # several stems on the same chain but no pair enclosing them
            return 'external'
        if not self.children or len(self.children) == 0:
            return 'hairpin'
        if len(self.children) == 1:
            return 'internal'
        if len(self.children) >= 2:
            return 'J' + str(len(self.children)+1)

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
        """
        This is the method that identifies the loop that corresponds to each node
        in the tree of nested basepairs.
        It is a Node method.
        It identifies the unpaired positions in each strand and stores the first
        index and the last index plus one for each strand.
        When the interior of a strand is empty, it stores a tuple like (14,14),
        because range(14,14) is empty.
        Make sure to avoid passing back an IL with two empty strands, because that is just two
        successive Watson-Crick pairs.
        But OK to pass back a J3 or J4 with all empty strands, because sometimes that happens.
        """

        unpaired = []
        left = self.left()
        found_unpaired_nts = False

        # print(self.children)  # just Node objects

        # sort the children by lowest index for J3, J4
        for child in sorted(self.children):
            right = child.value[0]
            # print('  Strand unpaired positions %s to %s, not including %s' % (left,right,right))

            # new: store the left and right, because in case they are equal,
            # then you can still tell where the flanking bases are supposed to be
            unpaired.append((left,right))

            if left < right:
                found_unpaired_nts = True
            # looped = range(left, right)
            # if looped:
            #     print('  Strand unpaired positions %s to %s appended' % (left,right))
            #     unpaired.append(looped)
            left = child.value[1] + 1
        # print('  Strand unpaired positions %s to %s, not including %s' % (left,self.value[1],self.value[1]))
        # last = range(left, self.value[1])
        # if last:
        #     print('  Strand unpaired positions %s to %s appended' % (left,self.value[1]))
            # unpaired.append(last)

        # new: store the left and right positions, in case the range is empty
        unpaired.append((left, self.value[1]))
        if left < self.value[1]:
            found_unpaired_nts = True

        if found_unpaired_nts or not len(self.children) == 1:
            # at least one strand has a non-empty interior,
            # or this is an HL, J3, J4, ... but not IL or external loop
            # print('basic.py: These are the children of node with left position %s' % unpaired[0][0])
            # for left, right in unpaired:
            #     print('  Strand unpaired positions %s to %s, not including %s' % (left,right,right))

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
        # print('basic gt: %s %s' % (self.value, other.value))
        # Python 3 does not allow comparison of None with integers
        # This craziness replicates the previous behavior
        if self.value is None:
            return False
        elif self.value[0] is None:
            return False
        elif self.value[1] is None:
            return False
        elif other.value[0] is None:
            return True
        elif other.value[1] is None:
            return True
        else:
            return self.value > other.value

    def __le__(self, other):
        return self.value <= other.value

    def __lt__(self, other):
        return self.value < other.value

    def __eq__(self, other):
        return isinstance(other, Node) and self.value == other.value and \
            self.children == other.children
