# This is an example of how to parse a dot-bracket string and get the loops.

from os import path
import sys

# Just mess with path a little so we can run this from anywhere.
here = path.abspath(path.join(path.dirname(__file__), '..'))
sys.path.append(here)

from rnastructure.secondary import dot_bracket as Dot

dot_string = "((((((..(((...)))...))))))"
sequence = "ccggccaacccuuugggcagggccgg"

parser = Dot.Parser(dot_string)

# Get all loop indices
print(parser.loops())
# => {'internal': [([6, 7], [17, 18, 19])], 'hairpin': [([11, 12, 13],)]}

# Get all loop indices plus flanking pairs
print(parser.loops(flanking=True))
# => {'internal': [([5, 6, 7, 8], [16, 17, 18, 19, 20])], 'hairpin': [([10, 11, 12, 13, 14],)]}

# Extact the sequences of all loops from a sequence. Note that the sequence and
# the structure must be the same length.
print(parser.parse(sequence))
# => {'internal': ['aa*cag'], 'hairpin': ['uuu']}

# Extract the sequences of all loops plus flanking pairs.
print(parser.parse(sequence, flanking=True))
# => {'internal': ['caac*gcagg'], 'hairpin': ['cuuug']}

# Other parsers are similar. The BPSeq and Connect parsers build with an
# iterable where each entry is a single line to parse. As an example an open
# file handle or a file which was slurped up and each line placed in a array.
# The methods to extract loops and indices are the same.
