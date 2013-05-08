""" This package contains a simple wrapper around the PDBx tools provided by
PDB.  It is intended to provide a simple and pythonic way to interface with
mmCIF data.
"""

import re
from itertools import islice

from pdbx.reader.PdbxParser import PdbxReader as Reader


class MissingBlockException(Exception):
    """This class is raised when trying to get a missing block of data.
    """
    pass


class MissingColumn(Exception):
    """This is raised when trying to get a missing column from a table.
    """
    pass


class CIF(object):
    """Top level container for all CIF related data. This assumes that each
    mmCIF file contains a single datablock. This doesn't have to be true but
    makes things easier.
    """

    def __init__(self, handle):
        reader = Reader(handle)
        self.data = []
        reader.read(self.data)
        self.data = self.data[0]

    def chain_polymer(self, requested):
        return [p for p in self.polymers() if p.chain == requested]

    def polymers(self):
        rows = []
        rows = list(self.pdbx_poly_seq_scheme.rows)
        chain = rows[0]['asym_id']
        for row in rows:
            monomer = row['auth_mon_id']
            if monomer == '?' or row['asym_id'] != chain:
                if rows:
                    yield Polymer(chain, rows)
                rows = []
                chain = row['asym_id']
            if monomer != '?':
                rows.append(row)
        if rows:
            yield Polymer(chain, rows)

    def polymer_sequences(self):
        for chain, polymer in self.polymers():
            sequence = [unit['mon_id'] for unit in polymer]
            yield chain, sequence

    def table(self, name):
        block_name = re.sub('^_', '', name)
        block = self.data.getObj(block_name)
        if not block:
            raise MissingBlockException("Unknown block " + name)
        return Table(block)

    def __getattr__(self, name):
        try:
            return self.table(name)
        except MissingBlockException:
            raise AttributeError("Unknown block " + name)


class Table(object):
    """Container for a single table in the data block. This provides some
    useful methods for accessing the data.
    """

    def __init__(self, block):
        self._block = block
        self.name = self._block.getName()
        self.columns = self._block.getItemNameList()
        self.columns = [re.sub('_.+\.', '', name) for name in self.columns]

    def row(self, number):
        """Get a row by index. Note that this may or may not be in the same
        order as they appear in the cif file, since cif files are not required
        to be ordered. The row will be a dict of the form { attribute: value }.
        Each attribute will have the name of the block stripped.
        """
        if number >= len(self):
            raise IndexError("Row index out of range.")

        return dict(zip(self.columns, self._block.getRow(number)))

    @property
    def rows(self):
        """Get a list of all rows"""
        for index in xrange(len(self)):
            yield self.row(index)

    def column(self, name):
        """Get a column by name"""
        if name not in self.columns:
            raise MissingColumn("Unknown column.")

        values = []
        for row in self.rows:
            values.append(row[name])
        return values

    def size(self):
        """Get a tuple of (rowCount, columnCount).
        """
        return (len(self), len(self.columns))

    def __getattr__(self, name):
        """Get the column with the given name.
        """
        try:
            return self.column(name)
        except MissingColumn:
            raise AttributeError("Unknown column")

    def __getitem__(self, index):
        """Get the row with the given index.
        """
        if isinstance(index, int):
            return self.row(index)
        if isinstance(index, str):
            try:
                return self.column(index)
            except MissingColumn:
                raise KeyError("Unknown column " + index)
        if isinstance(index, slice):
            # TODO: It would be nice to get another table back after slicing.
            iterator = islice(self.rows, index.start, index.stop, index.step)
            return list(iterator)
        raise TypeError("Unknown key type, should be str or int")

    def __len__(self):
        """Get the number of rows.
        """
        return self._block.getRowCount()


class TableSubset(Table):
    """This is used to represent a subset of a table. This will provide the
    normal Table methods but will use a list of rows.
    """
    def __init__(self, name, rows):
        self.name = 'pdbx_poly_seq_scheme'
        self._rows = rows
        self.columns = rows[0].keys()

    @property
    def rows(self):
        return self._rows

    def row(self, number):
        """Get a row by index. Note that this may or may not be in the same
        order as they appear in the cif file, since cif files are not required
        to be ordered. The row will be a dict of the form { attribute: value }.
        Each attribute will have the name of the block stripped.
        """
        return self.rows[number]

    def __len__(self):
        """Get the number of rows.
        """
        return len(self._rows)


class Polymer(TableSubset):
    """This represents a polymer in a mmCIF file. A polymer is a continuous set
    of resolved nucleotides from a single asymmetric unit.
    """
    def __init__(self, chain, rows):
        self.chain = chain
        super(Polymer, self).__init__('pdbx_poly_seq_scheme', rows)

    @property
    def sequence(self):
        return self.mon_id
