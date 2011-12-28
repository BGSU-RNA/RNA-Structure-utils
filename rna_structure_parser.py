from itertools import takewhile, ifilter, dropwhile


class RNAStructureParser:
    """A Simple class to parse 2D structures.
    This cannot handle pseudoknots.
    """

    def __init__(self, structure):
        stack = []
        pairs = [False] * len(structure)
        self.loops = { 'hairpins': [], 'internal': [] }

        for index, char in enumerate(structure):
            if char == '(':
                stack.append(index)
            elif char == ')':
                left = stack.pop()
                pairs[left] = tuple([index, '(', left])
                pairs[index] = tuple([left, ')', index])
            elif char == '.':
                pairs[index] = tuple([False, '.', index])
            else:
                raise "Unknown character"

        # Partition into each section
        helicies = []
        while self.__first_index(lambda p: p[1] == '(', pairs):
            helicies.append(self.partition(pairs))

        loops = []
        for helix in helicies:
            while helix:
                loop = self.__outer_loop(helix)
                if loop:
                    if helix:
                        self.loops['internal'].append(loop)
                    else:
                        self.loops['hairpins'].append(tuple(loop[0]))
                self.__strip_outer_helix(helix)

        # Get hairpin loops for each section

        # Link all internal loops

    def __outer(self, structure, index, check):
        outer = []
        while structure and check(structure[index]):
            outer.append(structure[index][2])
            del structure[index]
        return outer

    def __outer_loop(self, structure):
        left = self.__outer(structure, 0, lambda l: l[1] == '.')
        right = self.__outer(structure, -1, lambda l: l[1] == '.')
        if left or right:
            return (left, right)
        return False

    def __strip_outer_helix(self, structure):
        while structure and self.__pairs(structure[0], structure[-1]):
            del structure[0]
            del structure[-1]

    def __pairs(self, first, second):
        return first[0] == second[2] and second[2] == first[0]

    def __last_index(self, fn, structure):
        data = self.__first_index(fn, reversed(structure))
        if data:
            (i, f, c, l) = data
            return (len(structure) - i, f, c, l)
        return False

    def __first_index(self, fn, structure):
        enum = ifilter(lambda (i, p): fn(p), enumerate(structure))
        try:
            (i, (f, c, l)) = enum.next()
            return (i, f, c, l)
        except StopIteration:
            return False

    def partition(self, structure):
        (start, f, _, l) = self.__first_index(lambda p: p[1] == '(', structure)
        (stop, _, _, _) = self.__last_index(lambda p: p[2] == f, structure)
        parts = (structure[start:stop])
        del structure[start:stop]
        return parts

    def parse(self, sequence):
        ranges = {}
        for name, total in self.loops.iteritems():
            ranges[name] = []
            for positions in total:
                seq = ''.join(map(lambda v: sequence[v], positions))
                ranges[name].append(seq)
        return ranges


if __name__ == "__main__":
    parser = RNAStructureParser(".(((((.(((((((....................))))))))))..)).....(((((((((((((..............)))))))))))))....")
