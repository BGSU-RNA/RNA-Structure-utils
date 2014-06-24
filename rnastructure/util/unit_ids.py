import re

FRAGMENTS = ['pdb', 'model', 'chain', 'residue', 'number', 'atom_name',
             'alt_id', 'insertion_code', 'symmetry_operator']

SEPERATOR = '|'


class ImpossibleUnitIdException(Exception):

    """This is raised if we have a unit id that is totally invalid, for
    example, one that is empty or if the PDB is empty.
    """
    pass


class UnitIdGenerator(object):

    def __call__(self, obj, short=True):
        data = []
        for part in FRAGMENTS:
            if part in obj and obj[part] is not None:
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
        return SEPERATOR.join([d or '' for d in data])


class UnitIdParser(object):

    def __call__(self, unit_id):
        parts = unit_id.split(SEPERATOR)
        data = {name: None for name in FRAGMENTS}
        for index, part in enumerate(parts):
            name = FRAGMENTS[index]
            data[name] = part
        return data


def matlab_id_as_unit_id(pdb, matlab_id, **kwargs):
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


def as_matlab_id(unit_id):
    parser = UnitIdParser()
    data = parser(unit_id)
    residue = [data['residue'],
               str(data['number']),
               data['insertion_code'] or '']
    residue = ''.join(residue)
    return '%s:%s' % (data['chain'], residue)
