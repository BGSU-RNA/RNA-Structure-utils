import re

UNIT_ID = {
    'fragments': ['pdb', 'model', 'chain', 'residue', 'number', 'atom_name',
                  'alt_id', 'insertion_code', 'symmetry_operator'],
    'separator': '|',
    'required': ['pdb']
}

NUCLEOTIDE_ID = {
    'fragments': ['pdb', 'type', 'model', 'chain', 'number', 'residue',
                  'insertion_code'],
    'separator': '_',
    'required': ['pdb', 'type', 'model', 'chain', 'number', 'residue'],
}

MATLAB_ID = {
    'fragments': ['chain', 'residue', 'number', 'insertion_code'],
    'separator': (':', '', ''),
    'required': ['chain', 'residue', 'number'],
}


class ImpossibleNucleotideIdException(Exception):
    """This is raised if we have a nucleotide id that is totally invalid, for
    example, one that is empty or if the PDB is empty.
    """
    pass


class InvalidFragment(Exception):
    """This is raised if we are given a fragment name which is invalid.
    """
    pass


class MissingRequiredFragment(Exception):
    """This is raised when we have attempted to generate an id but we are
    missing a required fragment. For example, missing the chain in a matlab id,
    or pdb in a unit id.
    """
    pass


class InvalidUnitType(Exception):
    """This is raised when we are asked to build a generator with a bad name.
    """
    pass


class GenericIdGenerator(object):

    def __init__(self, config):
        self._config = config
        self._lookups = {}

    def insertion_code(self, data):
        if 'insertion_code' not in data:
            return None
        code = data['insertion_code']
        if code == '?' or not code or not code.strip():
            return None
        return code.strip()

    def fragments(self, obj, **kwargs):
        merged = {}
        for fragment in self._config['fragments']:
            if fragment in obj:
                merged[fragment] = obj[fragment]
            elif fragment in self._lookups:
                merged[fragment] = self._lookups[fragment](obj, **kwargs)
            else:
                merged[fragment] = kwargs.get(fragment, None)

        self.check_required(merged)
        data = []
        for fragment in self._config['fragments']:
            if fragment == 'insertion_code':
                data.append(self.insertion_code(merged))
            elif fragment not in merged:
                data.append(None)
            elif merged[fragment] is None:
                data.append(None)
            else:
                value = str(merged[fragment])
                if fragment == 'pdb':
                    value = value.upper()
                data.append(value)
        return data

    def check_required(self, data, **kwargs):
        for required in self._config['required']:
            if required not in data or data[required] is None:
                raise MissingRequiredFragment("Missing required %s" % required)
            if data[required] == '':
                raise MissingRequiredFragment("Given empty %s" % required)
        return True

    def lookup(self, name, func):
        self._lookups[name] = func


class UnitIdGenerator(GenericIdGenerator):

    def __init__(self):
        super(UnitIdGenerator, self).__init__(UNIT_ID)

    def __call__(self, obj, short=True, **kwargs):
        data = self.fragments(obj, **kwargs)

        # Check that both residue level entries are either set or not set
        if bool(data[3]) != bool(data[4]):
            raise MissingRequiredFragment("Must set both residue level ids")

        if short:

            # If no symmetry_operator or the default one, then strip it
            if data[-1] is None or data[-1] == '1_555':
                data = data[:-1]

            # Trim out as much as we can
            while data[-1] is None:
                data = data[:-1]

            # We should fail if we have chain but not model, or chain without
            # residue level stuff.
            for index, fragment in enumerate(data):
                if fragment is None and index < 3:
                    msg = 'Missing required %s' % UNIT_ID['fragments'][index]
                    raise MissingRequiredFragment(msg)

        # Sometimes data may be None, so we or it with the empty string to
        # get a string in all cases.
        return UNIT_ID['separator'].join([d or '' for d in data])


class UnitIdParser(object):

    def __call__(self, unit_id):
        parts = unit_id.split(UNIT_ID['separator'])
        data = dict((name, None) for name in UNIT_ID['fragments'])
        for index, part in enumerate(parts):
            name = UNIT_ID['fragments'][index]
            data[name] = part
        return data


class MatlabIdParser(object):

    def __call__(self, matlab_id, **kwargs):
        parts = matlab_id.split(':')
        number = parts[1][1:]
        ins = None

        if not re.match('\d', number[-1]):
            ins = number[-1]
            number = number[:-1]

        data = {
            'model': 1,
            'chain': parts[0],
            'residue': parts[1][0],
            'number': number,
            'insertion_code': ins
        }
        data.update(kwargs)
        return data


class MatlabIdGenerator(GenericIdGenerator):
    """This class is used to generate Matlab style ids for nucleotides.
    """

    def __init__(self):
        super(MatlabIdGenerator, self).__init__(MATLAB_ID)

    def __call__(self, data, **kwargs):
        raw = self.fragments(data, **kwargs)
        if raw[-1] is None:
            raw.pop(-1)
        return ':'.join([raw.pop(0), ''.join(raw)])


class NucleotideIdParser(object):
    def __call__(self, nt_id):
        parts = nt_id.split(NUCLEOTIDE_ID['separator'])
        data = dict(zip(NUCLEOTIDE_ID['fragments'], parts))
        if data['insertion_code'] == '':
            data['insertion_code'] = None
        return data


class NucleotideIdGenerator(GenericIdGenerator):
    def __init__(self):
        super(NucleotideIdGenerator, self).__init__(NUCLEOTIDE_ID)

    def __call__(self, data, **kwargs):
        fragments = self.fragments(data, **kwargs)
        if fragments[-1] is None:
            fragments[-1] = ''
        return NUCLEOTIDE_ID['separator'].join(fragments)


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


def generate_converter(input_type, output_type, **kwargs):
    parser = None
    if input_type == 'unit':
        parser = UnitIdParser()
    elif input_type == 'matlab':
        parser = MatlabIdParser()
    elif input_type == 'nucleotide' or input_type == 'old' \
            or input_type == 'nt':
        parser = NucleotideIdParser()
    else:
        raise InvalidUnitType("Input type %s is unknown" % input_type)

    writer = None
    if output_type == 'unit':
        writer = UnitIdGenerator()
    elif output_type == 'matlab':
        writer = MatlabIdGenerator()
    elif output_type == 'nucleotide' or output_type == 'old' \
            or input_type == 'nt':
        writer = NucleotideIdGenerator()
    else:
        raise InvalidUnitType("Output type %s is unknown" % output_type)

    return lambda string, **kwargs: writer(parser(string), **kwargs)
