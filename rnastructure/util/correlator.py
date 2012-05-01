from itertools import ifilter


def columns_to_index(sequence, names=None):
    size = len(sequence) - sequence.count('-')
    if names == None:
        names = range(size)

    if len(names) != size:
        raise ValueError("Must give as many names as non gap characters")

    correlations = {}
    current = 0
    for column, char in enumerate(sequence):
        if char != '-':
            correlations[column] = names[current]
            current += 1
    return correlations


def correlated_columns(reference, sequence):
    """
      Determine which columns are correlated between the two sequences.
    """
    indecies = range(len(reference))
    correlations = ifilter(lambda i: reference[i] != '-', indecies)
    correlations = ifilter(lambda i: sequence[i] != '-', correlations)
    return list(correlations)


def correlate_aligned_indices(reference, sequence, reference_names=None,
                              sequence_names=None):
    """
      Correlate the aligned indices in two sequences
    """
    ref_correlations = columns_to_index(reference, reference_names)
    seq_correlations = columns_to_index(sequence, sequence_names)
    correlated = correlated_columns(reference, sequence)
    correlations = {}
    for column in correlated:
        ref_column = ref_correlations[column]
        seq_column = seq_correlations[column]
        correlations[ref_column] = seq_column
    return correlations
