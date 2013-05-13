"""
A simple container for sequences. Basically just a way to assign names to
positions in the sequence.
"""


class Sequence(str):
    def __init__(self, sequence, names=None):
        if names is None:
            count = len(sequence) - sequence.count('-')
            names = range(count)
        self.names = names
        super(Sequence, self).__init__(sequence)
