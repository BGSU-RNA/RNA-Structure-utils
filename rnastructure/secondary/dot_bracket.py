from rnastructure.secondary.basic import Basic


class Parser(Basic):
    """A class to parse 2D structures in dot-bracket format.

    The dot bracket format must be composed of <, (, [, {, ., :, }, ], ), >
    only. Paired bases are represented by (, <, >, and ). Unpaired bases are
    represented by ., and : while pseudonoted characters are represented by {,
    [, ], and }.

    """

    def __init__(self, structure):
        """Construct a new Parser object with the given structure.

        The structure should be a string in dot bracket notation with possible
        pseudoknots.

        """
        pairs = self.pairs(structure)
        super(Parser, self).__init__(pairs)

    def pairs(self, structure):
        """Compute which bases are paired in the 2D structure.
        """
        helix_stack = []
        knot_stack = []
        pairs = [None] * len(structure)
        for index, char in enumerate(structure):
            if char == '(' or char == '<':
                helix_stack.append(index)
            elif char == ')' or char == '>':
                left = helix_stack.pop()
                pairs[left] = index
                pairs[index] = left
            elif char == '.' or char == ':':
                pass
            elif char == '{' or char == '[':
                knot_stack.append(index)
            elif char == '}' or char == ']':
                left = knot_stack.pop()
                pairs[left] = index
                pairs[index] = left
            else:
                raise ValueError("Unknown character: '%s'" % char)
        return pairs
