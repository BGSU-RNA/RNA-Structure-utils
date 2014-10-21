#!/usr/bin/env python

from os import path
import sys

import argparse

# Just mess with path a little so we can run this from anywhere.
here = path.abspath(path.join(path.dirname(__file__), '..'))
sys.path.insert(0, here)

from rnastructure.tertiary.cif import CIF


def main(cif, symm, model, chain):
    with open(cif, 'rb') as raw:
        structure = CIF(raw)

    chain = structure.chain(symm, model, chain)
    for seq, seq_id, unit_id in chain.experimental_sequence_mapping():
        print('%s\t%s' % (seq_id, unit_id))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("cif", help="Cif file to read")
    parser.add_argument("--chain", default="A",
                        help="Chain to use")
    parser.add_argument("--model", default=1,
                        help="Model to use")
    parser.add_argument('--symmetry', default='1_555',
                        help="Symmetry operator to use")
    args = parser.parse_args()

    main(args.cif, args.symmetry, args.model, args.chain)
