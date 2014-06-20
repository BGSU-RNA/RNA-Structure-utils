""" This package contains a simple wrapper around the PDBx tools provided by
PDB.  It is intended to provide a simple and pythonic way to interface with
mmCIF data.
"""

import re
import functools as func
import collections as coll
import itertools as it

from pdbx.reader.PdbxParser import PdbxReader as Reader


class MissingBlockException(Exception):

    """This class is raised when trying to get a missing block of data.
    """
    pass


class MissingColumn(Exception):

    """This is raised when trying to get a missing column from a table.
    """
    pass


class ComplexOperatorException(Exception):

    """This is raised when we come across complex operators that we cannot
    easily deal with. These tend to show up in viral structures and not things
    we deal with currently.
    """
    pass


class UnusableUnobservedTable(Exception):
    pass


class ImpossibleUnitIdException(Exception):

    """This is raised if we have a unit id that is totally invalid, for
    example, one that is empty or if the PDB is empty.
    """
    pass


def atom_sorter(atom):
    return (atom['label_asym_id'], int(atom['pdbx_PDB_model_num']),
            atom['auth_asym_id'], int(atom['auth_seq_id']))


def filter_using(iterable, key, target):
    return it.ifilter(lambda o: o[key] == target, iterable)


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
        self.name = self.data.getName()
        self._assemblies = self.__load_assemblies__()
        self._entities = self.__load_entities__()

    def __load_assemblies__(self):
        operators = {op['id']: op for op in self.pdbx_struct_oper_list}
        assemblies = coll.defaultdict(list)
        for assembly in self.pdbx_struct_assembly_gen:
            operator = assembly['oper_expression']
            if operator not in operators:
                raise ComplexOperatorException()
            for asym_ids in assembly['asym_id_list'].split(','):
                for asym_id in asym_ids:
                    assemblies[asym_id].append(operators[operator])
        return assemblies

    def __load_entities__(self):
        entities = {}
        for entity in self.entity:
            entities[entity['id']] = entity
        return entities

    def symmetry_operators(self, **kwargs):
        atoms = sorted(self.atom_site.rows, key=atom_sorter)
        for operator in self.pdbx_struct_oper_list:
            fn = lambda a: operator in self.operators(a['label_asym_id'])
            yield Symmetry(self, operator, it.ifilter(fn, atoms), **kwargs)

    def symmetry_operator(self, name, **kwargs):
        operator = None
        for row in self.pdbx_struct_oper_list:
            if row['name'] == name:
                operator = row
                break

        if not operator:
            return None

        atoms = sorted(self.atom_site.rows, key=atom_sorter)
        fn = lambda a: operator in self.operators(a['label_asym_id'])
        op = Symmetry(self, operator, it.ifilter(fn, atoms), **kwargs)

        if not op:
            return None
        return op

    def models(self):
        for op in self.symmetry_operators():
            for model in op.models():
                yield model

    def model(self, operator, model_number):
        op = self.symmetry_operator(operator)
        if not op:
            return None
        return op.model(model_number)

    def chains(self):
        for model in self.models():
            for chain in model.chains():
                yield chain

    def chain(self, operator, model_number, chain_id):
        model = self.model(operator, model_number)
        if not model:
            return None
        return model.chain(chain_id)

    def polymers(self):
        for chain in self.chains():
            for polymer in chain.polymers():
                yield polymer

    def table(self, name):
        return Table(self, self.__block__(name))

    def operators(self, asym_id):
        return self._assemblies[asym_id]

    def is_water(self, entity_id):
        return self._entities[entity_id]['type'] == 'water'

    def is_polymeric(self, entity_id):
        return self._entities[entity_id]['type'] == 'polymer'

    def is_polymeric_atom(self, atom):
        return self.is_polymeric(atom['label_entity_id'])

    def __block__(self, name):
        block_name = re.sub('^_', '', name)
        block = self.data.getObj(block_name)
        if not block:
            raise MissingBlockException("Unknown block " + name)
        return block

    def __getattr__(self, name):
        try:
            return self.table(name)
        except MissingBlockException:
            raise AttributeError("Unknown block " + name)


class Table(object):

    """Container for a single table in the data block. This provides some
    useful methods for accessing the data.
    """

    def __init__(self, cif, block, rows=None):
        self._cif = cif
        self.block = block
        self.rows = rows

        self.columns = self.block.getItemNameList()
        self.columns = [re.sub('_.+\.', '', name) for name in self.columns]

        if self.rows is None:
            length = self.block.getRowCount()
            self.rows = [self.__row__(index) for index in xrange(length)]

    def column(self, name):
        """Get a column by name"""
        if name not in self.columns:
            raise MissingColumn("Unknown column")

        values = []
        for row in self.rows:
            values.append(row[name])
        return values

    def size(self):
        """Get a tuple of (rowCount, columnCount).
        """
        return (len(self), len(self.columns))

    def __row__(self, number):
        """Get a row by index. Note that this may or may not be in the same
        order as they appear in the cif file, since cif files are not required
        to be ordered. The row will be a dict of the form { attribute: value }.
        Each attribute will have the name of the block stripped.
        """
        return dict(zip(self.columns, self.block.getRow(number)))

    def __getattr__(self, name):
        """Get the column with the given name.
        """
        try:
            return self.column(name)
        except MissingColumn:
            raise AttributeError("Unknown column: %s" % name)

    def __getitem__(self, index):
        if isinstance(index, str):
            try:
                return self.column(index)
            except MissingColumn:
                raise KeyError("Unknown column: %s" % index)

        if isinstance(index, int):
            return self.rows[index]

        if isinstance(index, slice):
            return Table(self._cif, self.block, rows=self.rows[index])

        raise TypeError("Unknown key type, should be str, int or slice")

    def __len__(self):
        """Get the number of rows.
        """
        return len(self.rows)


class GenericMapping(coll.Mapping):
    _properties = {}

    def inherit(self, obj, **kwargs):
        for key, value in obj.items():
            self._properties[key] = value
        self._properties.update(kwargs)

    def __getitem__(self, key):
        return self._properties[key]

    def __iter__(self):
        return iter(self._properties)

    def __len__(self):
        return len(self._properties)


class ResidueContainer(object):

    def __init__(self, cif, atoms, unobs=None, nonpolymers=False):
        self._cif = cif

        if not nonpolymers:
            atoms = it.ifilter(cif.is_polymeric_atom, atoms)

        if not isinstance(atoms, list):
            atoms = list(atoms)

        self._atoms = atoms
        self._unobs = unobs

    @property
    def unobs(self):
        if self._unobs is None:
            unobs = self._cif.pdbx_unobs_or_zero_occ_residues.rows

            if 'chain' in self:
                unobs = filter_using(unobs, 'auth_asym_id', self['chain'])

            if 'model' in self:
                unobs = filter_using(unobs, 'PDB_model_num', self['model'])

            unobs = it.ifilter(lambda u: u['polymer_flag'] == 'Y', unobs)
            self._unobs = sorted(unobs, key=lambda u: int(u['auth_seq_id']))

        return self._unobs

    def atoms(self):
        for atom in self._atoms:
            yield atom

    def residues(self):
        for residue, atoms in self.__grouped__():
            yield Residue(self._cif, self['symmetry_operator'], list(atoms))

    def residue(self, target):
        for index, residue in enumerate(self.residues()):
            if index == target:
                return residue
        raise IndexError()

    def __grouped__(self):
        fn = lambda r: (r['auth_seq_id'], r['pdbx_PDB_ins_code'])
        return it.groupby(self.atoms(), fn)

    def __bool__(self):
        return bool(self._atoms)

    __nonzero__ = __bool__

    def __len__(self):
        count = 0
        for _, _ in self.__grouped__():
            count += 1
        return count


class UnitIdGenerator(object):

    part_ordering = ['pdb', 'model', 'chain', 'residue', 'number', 'atom_name',
                     'alt_id', 'insertion_code', 'symmetry_operator']
    seperator = '|'

    def __call__(self, obj, short=True):
        data = []
        for part in self.part_ordering:
            if part in obj:
                data.append(str(obj[part]))
            else:
                data.append(None)

        # Convert ? to None for insertion code, since that means not present in
        # cif files
        if data[7] == '?':
            data[7] = None

        if data[0] is None:
            raise ImpossibleUnitIdException("Can't make Unit id without PDB")

        # Check that both residue level entries are either set or not set
        if bool(data[3]) != bool(data[4]):
            raise ImpossibleUnitIdException("Must set both residue level ids")

        if short:
            # If no symmetry_operator or the default one, then strip it
            if data[-1] is None or data[-1] == '1_555':
                data = data[:-1]

            # Trim out as much as we can
            while data[-1] is None:
                data = data[:-1]

        # Sometimes data may be None, so we or it with the empty string to
        # get a string in all cases.
        return self.seperator.join([d or '' for d in data])


class Symmetry(ResidueContainer, GenericMapping):

    def __init__(self, cif, operator, atoms):
        super(Symmetry, self).__init__(cif, atoms)
        self.inherit({'pdb': cif.name, 'symmetry_operator': operator['name']})

    def model(self, number, **kwargs):
        num = str(number)
        fn = lambda r: r['pdbx_PDB_model_num'] == num
        atoms = it.ifilter(fn, self.atoms())
        model = Model(self._cif, num, self, atoms, **kwargs)
        if not model:
            return None
        return model

    def models(self, **kwargs):
        fn = lambda r: r['pdbx_PDB_model_num']
        for model_number, atoms in it.groupby(self.atoms(), fn):
            model = Model(self._cif, model_number, self, atoms, **kwargs)
            if model:
                yield model


class Model(ResidueContainer, GenericMapping):

    def __init__(self, cif, model_number, operator, atoms):
        super(Model, self).__init__(cif, atoms)
        self.unit_id = func.partial(UnitIdGenerator(), self)
        self.inherit(operator, model=model_number)

    def chain(self, chain_id, **kwargs):
        fn = lambda r: r['auth_asym_id'] == chain_id
        atoms = it.ifilter(fn, self.atoms())
        chain = Chain(self._cif, chain_id, self, atoms, **kwargs)
        if not chain:
            return None
        return chain

    def chains(self, **kwargs):
        grouped = it.groupby(self.atoms(), lambda r: r['auth_asym_id'])
        for chain_id, atoms in grouped:
            chain = Chain(self._cif, chain_id, self, atoms, **kwargs)
            if chain:
                yield chain


class Chain(ResidueContainer, GenericMapping):

    def __init__(self, cif, chain_id, model, atoms, **kwargs):
        super(Chain, self).__init__(cif, atoms, **kwargs)
        self.inherit(model, chain=chain_id)
        self.unit_id = func.partial(UnitIdGenerator(), self)

    def polymers(self):
        """Creates an iterator over each part of the chain which is a polymer.
        That means it is connected and marked as a polymer in the file. Any
        monomers that are unoboserved or have no occupancy cause a chain break.
        """

        def create_filter(endpoint):
            number = int(endpoint['auth_seq_id'])

            def func(atom):
                atom_number = int(atom['auth_seq_id'])
                if atom_number < number:
                    return True
                if atom_number == number:
                    return endpoint['PDB_ins_code'] == '?' or \
                        atom['PDB_ins_code'] < endpoint['PDB_ins_code']
                return False

            return func

        atoms = self.atoms()
        for index, endpoint in enumerate(self.unobs):
            if 'auth_seq_id' not in endpoint:
                raise UnusableUnobservedTable()

            polymer = list(it.takewhile(create_filter(endpoint), atoms))
            atoms = it.dropwhile(create_filter(endpoint), atoms)

            unobs = self.unobs[index + 1:]
            chain = Chain(self._cif, self['chain'], self, polymer, unobs=unobs)
            if chain:
                yield chain

        chain = Chain(self._cif, self['chain'], self, atoms, unobs=[])
        if chain:
            yield chain

    def polymer(self, target):
        for index, polymer in enumerate(self.polymers()):
            if index == target:
                return polymer
        raise IndexError()

    def has_breaks(self):
        first_seq_id = int(self._atoms[0]['auth_seq_id'])
        last_seq_id = int(self._atoms[-1]['auth_seq_id'])
        first_unobs = int(self.unobs[0]['auth_seq_id'])
        last_unobs = int(self.unobs[-1]['auth_seq_id'])
        return last_seq_id > first_unobs and first_seq_id < last_unobs

    @property
    def sequence(self):
        """The sequence of this chain as a array of three letter codes.
        """
        return [residue['residue'] for residue in self.residues()]


class Residue(GenericMapping):

    def __init__(self, cif, symmetry_operator, atoms):
        super(Residue, self).__init__()
        self.inherit({
            'pdb': cif.name,
            'model': atoms[0]['pdbx_PDB_model_num'],
            'chain': atoms[0]['auth_asym_id'],
            'residue': atoms[0]['label_comp_id'],
            'number': atoms[0]['auth_seq_id'],
            'insertion_code': atoms[0]['pdbx_PDB_ins_code'],
            'symmetry_operator': symmetry_operator,
            'residue': atoms[0]['auth_comp_id']
        })
        print(self._properties)
        self._cif = cif
        self._atoms = atoms
        self.unit_id = func.partial(UnitIdGenerator(), self)

    def atoms(self):
        for atom in self._atoms:
            yield atom

    def __len__(self):
        return len(self._atoms)
