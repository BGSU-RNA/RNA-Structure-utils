"""
This module has useful things for processing R3D Align output
"""

import csv

from rnastructure.util.unit_ids import matlab_id_as_unit_id as as_nt_id


class UnclearMappingException(Exception):
    """This is raised when we are not sure if how to correlate nucleotides.
    This occurs when we detect a two pairs where A maps to B as well as A
    mapping to C.
    """
    pass


HEADER = ['nt1', 'bp1', 'nt2', 'nt3', 'bp2', 'nt4', 'discrepancy', 'EMPTY']


def bp2nt(raw):
    raw_header = raw.readline().split(',')
    pdb1 = raw_header[1]
    pdb2 = raw_header[4]
    reader = csv.DictReader(raw, fieldnames=HEADER)
    data = []
    seen = {}
    for row in reader:
        if row['nt1'] == '---' or row['nt3'] == '---':
            continue
        nt1 = as_nt_id(pdb1, row['nt1'])
        nt2 = as_nt_id(pdb2, row['nt3'])
        if nt1 in seen and seen[nt1] != nt2:
            raise UnclearMappingException()
        seen[nt1] = nt2
        seen[nt2] = nt1
        data.append((nt1, nt2))
    return data
