import os
import sys

# Just mess with path a little so we can run this from anywhere.
here = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(here)

from rnastructure.secondary.pseudoknot import RemovePseudoknots
from rnastructure.secondary.dot_bracket import Parser as DotParser
from rnastructure.secondary.dot_bracket import Writer as DotWriter

rfam_pseudoknots = '<<<<<..AA..>>>>>aa'
rfam_sequence    = 'CCCCCAAGGUUGGGGGCC'

# Create a dot-bracket parser that can handle rfam's format.
parser = DotParser(rfam_pseudoknots, dialect='rfam')
parser.sequence = rfam_sequence

print("""
Note that RemovePseudoknots needs to have an enviroment variable DATAPATH set.
This must lead to the data tables provided with RNAstructure. Below is where
mine are, change as needed.
os.putenv('DATAPATH', '/Users/blake/lab/programs/RNAstructure/data_tables')
""")

# Remove all pseudoknots and return a new parser for the new structure.
remover = RemovePseudoknots()
stripped = remover(parser)

# Show the indices
print(stripped.indices())
# => {'external': [([16, 17],)], 'hairpin': [([5, 6, 7, 8, 9, 10],)]}

# Create a new dot-bracket string without pseudoknots - This does not support
# dialects at the moment, maybe later.
writer = DotWriter()
print(writer.format(stripped))
# => (((((......)))))..
