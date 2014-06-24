"""
This module has useful things for processing R3D Align output
"""

import re
import csv

from rnastructure.tertiary.cif import UnitIdGenerator


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


def as_nt_id(pdb, matlab_id, **kwargs):
    """Convert one of the matlab ids to something like a unit id. It is likely
    to be a correct unit id, however by default we are assuming that the
    symmetry operator is 1_555 and the model is 1. This can be changed with
    keyword arguments.
    """

    generator = UnitIdGenerator()
    parts = matlab_id.split(':')
    number = parts[1][1:]
    ins = None
    if not re.match('\d', number[-1]):
        ins = number[-1]
        number = number[:-1]

    data = {
        'pdb': pdb,
        'model': 1,
        'chain': parts[0],
        'residue': parts[1][0],
        'number': number,
        'insertion_code': ins
    }
    data.update(kwargs)
    return generator(data)
