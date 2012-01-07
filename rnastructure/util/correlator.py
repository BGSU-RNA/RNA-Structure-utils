class Correlator(object):

    def __init__(self, reference, names=None):
        if names == None:
            names = range(len(reference))

        self._names = names
        self._reference = reference

    def correlate(sequence, names=None):
        
