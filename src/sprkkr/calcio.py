# coding: utf-8
# vim: set et sw=4 ts=4 sts=4 nu ai cc=80 fdm=indent mouse=a:

"""
Module calcio
=============

"""
import os
import numpy as np
from datetime import datetime
from collections import OrderedDict
from ast import literal_eval

from ase import Atom
from ase.units import Bohr
from ase.spacegroup import get_spacegroup

from .misc import LOGGER, get_occupancy, AttrDict
from .data import get_sprkkr_input

_tr = str.maketrans('{}', '[]')

def _parse_value(val):
    aux = val.strip().translate(_tr)
    try:
        out = literal_eval(aux)

    except (ValueError, SyntaxError):
        aux2 = aux.split()
        if len(aux2) == 2: # Treat value with a unit.
            try:
                out = literal_eval(aux2[0])

            except (ValueError, SyntaxError):
                out = aux

        else:
            out = aux

    return out

def load_parameters(lines):
    """
    Load parameters.

    Notes
    -----
    The following assumptions should hold:
    - The section name is the first word on a line.
    - The sections are separated by empty lines or comments.
    - Each entry is on a separate line, with the exception of the first entry
      following the section name.

    Returns
    -------
    pars : dict
        The parameters as {section_name : section_dict}. When a parameter takes
        no value, the returned value is '@defined'.
    """
#    with open(filename) as fd:
#        lines = fd.readlines()

    # Replace comments.
    lines = [line.strip() if not line.startswith('#') else ''
             for line in lines]
    # Split sections.
    sections = '\n'.join(lines).split('\n\n')

    # Parse sections.
    pars = OrderedDict()
    for section in sections:
        section = section.strip()
        if not len(section): continue

        # Get section name.
        items = section.split('\n')
        name = items[0].split()[0]
        items[0] = items[0].replace(name, '')
        name = name.lower()

        # Parse section items.
        vals = OrderedDict()
        for item in items:
            aux = item.split('=')
            key = aux[0].strip()
            if len(aux) == 2: # item is key = value
                val = _parse_value(aux[1])

            else: # item is flag
                val = '@defined'

            vals[key] = val

        pars[name] = vals

    return pars

def make_sections(pars):
    sections = OrderedDict()
    for name, vals in pars.items():
        section = Section(name)

        for key, item in vals.items():
            dim = 0
            if isinstance(item, int):
                fmt = '{:d}'

            elif isinstance(item, float):
                fmt = '{:f}'

            elif isinstance(item, list):
                fmt = '{}'
                dim = len(item)

            else:
                fmt = '{:s}'

            var = Variable(key, default=item, fmt=fmt, dim=dim,
                           always_show=True)
            section.add_variables(var)

        sections[name] = section

    return sections

def compare_sections(s1, s2):
    is_eq = {name : (section == s2[name]) for name, section in s1.items()}
    return is_eq

class Section(object):
    def __init__(self, name):
        self.name = name
        self.variables = []
        self.switches = []

    def add_variables(self, *variables):
        # dynamically add an attribute
        for variable in variables:
            self.variables.append(variable)

    def set(self, **kwargs):
        names = [_.name.upper() for _ in self.variables]
        for kw, val in kwargs.items():
            if kw.upper() not in names:
                self.set_custom(**{kw:val})
            else:
                for v in self.variables:
                    if kw.upper() == v.name.upper():
                        v.value = val

    def set_custom(self, **kwargs):
        LOGGER.info(f"Adding custom variables to the \"{self.name}\" section")
        for kw, val in kwargs.items():
            v = Variable(kw)
            v.value = val
            self.add_variables(v)

    def __str__(self):
        # for pretty printing, get the indentation
        indent = len(self.name)
        s = self.name.upper()
        for i, variable in enumerate(self.variables):
            if i > 0:
                s += " " * indent
            s += f" {variable}\n"
        return s

    def __eq__(self, other):
        if len(self.variables) != len(other.variables):
            return False

        if len(self.switches) != len(other.switches):
            return False

        for ii, var in enumerate(self.variables):
            if var != other.variables[ii]:
                return False

        return True

class Variable(object):
    def __init__(self, name, default=None, fmt="{:s}", dim=0,
                 always_show=True):
        self.name = name
        self.fmt = fmt
        self.dim = dim
        self.default = default
        self.value = self.default
        self.always_show = always_show

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, val):
        self._value = val

    def __str__(self):
        if self.value is None:
            raise NameError(f"You must set a value for \"{self.name}\"")
        if self.value == self.default and not(self.always_show):
            return ""

        s   = val = ""
        key = self.name.upper()
        if self.dim > 0:
            try:
                val += ",".join([self.fmt.format(_) for _ in self.value])
            except:
                val += self.fmt.format(self.value)
            finally:
                s = key + "={" + val + "}"
        else:
            if self.value =='@defined':
               s = key
            else:
               s = key + "=" + self.fmt.format(self.value)
        return s

    def __eq__(self, other):
        is_eq = ((self.name == other.name)
                 and (self.value == other.value)
                 and (self.dim == other.dim))
        return is_eq

class InputFile(object):
    def __init__(self, filename="kkr.inp", task='scf',directory='./'):
        self.filename = filename
        self.directory=directory
        self.task=task
        self.set_defaults()

    def set_defaults(self):
        lines=get_sprkkr_input(self.task)
        if lines == 'None': # Empty sections - to remove?
            control_section = Section("CONTROL")
            sites_section = Section("SITES")
            tau_section = Section("TAU")
            scf_section=Section("SCF")

            self.sections = OrderedDict({
                "control": control_section,
                "tau"    : tau_section,
                "scf"    : scf_section,
                "sites"  : sites_section,
            })
        else:
            lines=get_sprkkr_input(self.task)
            pars = load_parameters(lines)
            self.sections = make_sections(pars)
        for key, section in self.sections.items():
            setattr(self, key + "_section", section)

    def set(self, section_name, **variables):
        if not(hasattr(self, section_name + "_section")):
            LOGGER.info(f"Adding custom section \"{section_name}\"")
            setattr(self, section_name + "_section", Section(section_name))

        section = getattr(self, section_name + "_section")
        section.set(**variables)
        self.sections["section_name"] = section

    def write(self):
        inp=os.path.join(self.directory, self.filename)
        with open(inp, "w") as fd:
            fd.write(str(self))

    def __str__(self):
        now = datetime.now()
        s  = "#"*80 + "\n"
        s += "#  SPR-KKR input file   {}\n".format(self.filename)
        s += "#  created by sprkkr Python package on {}\n".format(now.ctime())
        s += "#"*80 + "\n"
        s += "\n"

        for key, section in self.sections.items():
            s += str(section)
            s += "\n"

        return s

def _skip_lines(fd, num):
    for ii in range(num):
        line = next(fd)
    return line

def _skip_lines_to(fd, key):
    while 1:
        try:
            line = next(fd)

        except StopIteration:
            return ''

        if key in _dec(line):
            return line

def _dec(val):
    if isinstance(val, bytes):
        return val.decode('utf-8')

    else:
        return val

def _next_key_val(fd, dtype=bytes):
    aux = next(fd).split(maxsplit=1)
    return _dec(aux[0].lower()), dtype(aux[1].strip())

def _gen_next_key_val(fd, dtype=bytes):
    for line in fd:
        aux = line.split(maxsplit=1)
        if len(aux) == 2:
            yield _dec(aux[0].lower()), dtype(aux[1].strip())

        else:
            return

def _convert_array(x, dtype):
    return np.array([dtype(ii) for ii in _dec(x).split()])

_scf_info_dtypes = {
    int
    : {'ibzint', 'nktab', 'scf-iter'},
    float
    : {'scf-mix', 'scf-tol', 'rmsavv', 'rmsavb', 'ef', 'vmtz'},
    lambda x: _convert_array(x, int)
    : {'ne'},
}

_lattice_dtypes = {
    lambda x: [fun(x.split()[ii]) for ii, fun in enumerate([int] + 4 * [str])]
    : {'bravais'},
    float
    : {'alat'},
    lambda x: _convert_array(x, float)
    : {'a(1)', 'a(2)', 'a(3)'},
}

_sites_dtypes = {
    lambda x: _convert_array(x, float)
    : {'basscale'},
}

_reference_dtypes = {
    int
    : {'nref'},
}

_magdir_dtypes = {
    int
    : {'kmrot'},
    lambda x: _convert_array(x, float)
    : {'qmvec'},
}

def _convert(val, key, dtypes):
    for dtype, keys in dtypes.items():
        if key in keys:
            return dtype(val)

    else:
        return _dec(val)

def read_potential(filename):
    class PotFromFile(AttrDict):
        pass

    out = PotFromFile(filename=filename)

    with open(filename, 'rb') as fd:
        _skip_lines(fd, 1)
        setattr(out, *_next_key_val(fd))
        assert(out.get('header') is not None)

        _skip_lines(fd, 1)
        for _ in range(4):
            setattr(out, *_next_key_val(fd))
        assert(out.get('title') is not None)

        line = _skip_lines(fd, 2)
        assert('GLOBAL SYSTEM PARAMETER' in _dec(line))
        for _ in range(5):
            setattr(out, *_next_key_val(fd, dtype=int))

        line = _skip_lines(fd, 2)
        assert('SCF-INFO' in _dec(line))
        for key, val in _gen_next_key_val(fd):
            out[key] = _convert(val, key, _scf_info_dtypes)

        line = _skip_lines(fd, 1)
        assert('LATTICE' in _dec(line))
        for key, val in _gen_next_key_val(fd):
            out[key] = _convert(val, key, _lattice_dtypes)

        line = _skip_lines(fd, 1)
        assert('SITES' in _dec(line))
        for key, val in _gen_next_key_val(fd):
            if key == 'iq':
                out.qbas = np.loadtxt(fd, usecols=(1, 2, 3), max_rows=out.nq,
                                      ndmin=2)
                break
            out[key] = _convert(val, key, _sites_dtypes)

        line = _skip_lines(fd, 2)
        assert('OCCUPATION' in _dec(line))
        for key, val in _gen_next_key_val(fd):
            if key == 'iq':
                out.occupation = []
                for iq in range(out.nq):
                    line = next(fd).split()
                    aux1 = list(map(int, line[1:4]))
                    aux2 = [fun(line[4+ii])
                            for ii, fun in enumerate([int, float] * aux1[2])]
                    out.occupation.append(aux1 + [aux2])
                break

        line = _skip_lines(fd, 2)
        assert('REFERENCE SYSTEM' in _dec(line))
        for key, val in _gen_next_key_val(fd):
            if key == 'iref':
                out.ref = np.loadtxt(fd, usecols=(1, 2), max_rows=out.nref,
                                     ndmin=2)
                break
            out[key] = _convert(val, key, _reference_dtypes)

        line = _skip_lines(fd, 2)
        assert('HOST MADELUNG POTENTIAL' in _dec(line))
        line = _skip_lines_to(fd, 'NLMTOP-POT').split()
        out.nlmtop_pot_num = int(line[1])
        out.nlmtop_pot = np.loadtxt(
            fd, dtype=np.dtype('i4, i4, f8'),
            max_rows=out.nq * out.nlmtop_pot_num,
            ndmin=1,
        )

        line = _skip_lines(fd, 2)
        assert('CHARGE MOMENTS' in _dec(line))
        line = _skip_lines_to(fd, 'NLMTOP-CHR').split()
        out.nlmtop_chr_num = int(line[1])
        out.nlmtop_chr = np.loadtxt(
            fd, dtype=np.dtype('i4, i4, f8'),
            max_rows=out.nq * out.nlmtop_chr_num,
            ndmin=1,
        )

        line = _skip_lines(fd, 2)
        assert('MAGNETISATION DIRECTION' in _dec(line))
        # Ignores M*_T entries.
        for key, val in _gen_next_key_val(fd):
            if key == 'iq':
                aux = np.loadtxt(fd, max_rows=out.nq, ndmin=2)
                out.mag_dir = aux[:, 1:]
                break
            out[key] = _convert(val, key, _magdir_dtypes)

        _skip_lines_to(fd, 'MESH INFORMATION')
        setattr(out, *_next_key_val(fd))
        _skip_lines(fd, 1)
        out.mesh = np.loadtxt(
            fd, dtype=np.dtype([('im', 'i4'), ('r_1', 'f8'), ('dx', 'f8'),
                                ('jrmt', 'i4'), ('rmt', 'f8'), ('jrws', 'i4'),
                                ('rws', 'f8')]),
            max_rows=out.nm,
            ndmin=1,
        )

        line = _skip_lines(fd, 2)
        assert('TYPES' in _dec(line))
        line = _skip_lines(fd, 1).split()
        out.types = np.loadtxt(
            fd, dtype=np.dtype([('it', 'i4'), ('txt_t', 'U10'), ('zt', 'i4'),
                                ('ncort', 'i4'), ('nvalt', 'i4'),
                                ('nsemcorshlt', 'i4')]),
            max_rows=out.nt,
            ndmin=1,
        )

        line = _skip_lines(fd, 2)
        assert('POTENTIAL' in _dec(line))
        out.potentials_vt = []
        out.potentials_bt = []
        for ii in range(out.nt):
            line = _skip_lines(fd, 1).split()
            assert((_dec(line[0]) == 'TYPE') and
                   (int(line[1]) == out.types[ii]['it']))
            val = np.fromfile(fd, count=out.mesh[0]['jrws'], sep=' ',
                              dtype=np.float64)
            out.potentials_vt.append(val)
            val = np.fromfile(fd, count=out.mesh[0]['jrws'], sep=' ',
                              dtype=np.float64)
            out.potentials_bt.append(val)
            _skip_lines_to(fd, '========')

        line = _skip_lines(fd, 2)
        assert('CHARGE' in _dec(line))
        out.charges = []
        out.charges_bt = []
        for ii in range(out.nt):
            line = _skip_lines(fd, 1).split()
            assert((_dec(line[0]) == 'TYPE') and
                   (int(line[1]) == out.types[ii]['it']))
            val = np.fromfile(fd, count=out.mesh[0]['jrws'], sep=' ',
                              dtype=np.float64)
            out.charges.append(val)
            val = np.fromfile(fd, count=out.mesh[0]['jrws'], sep=' ',
                              dtype=np.float64)
            out.charges_bt.append(val)
            _skip_lines_to(fd, '========')

        line = _dec(_skip_lines(fd, 2)).split()
        assert('MOMENTS' in line[0])
        num = len(line[1:])
        out.moments_names = line[1:]
        out.moments = []
        for ii in range(out.nt):
            line = _skip_lines(fd, 1).split()
            assert((_dec(line[0]) == 'TYPE') and
                   (int(line[1]) == out.types[ii]['it']))
            val = np.fromfile(fd, count=num, sep=' ', dtype=np.float64)
            out.moments.append(val)
            _skip_lines_to(fd, '========')

    return out

class LatticeData(object):
    sym = OrderedDict({
            'aP': 'C_i',
            'mP': 'C_2h', 'mS': 'C_2h',
            'oP': 'D_2h', 'oS': 'D_2h', 'oI': 'D_2h', 'oF': 'D_2h',
            'tP': 'D_4h', 'tI': 'D_4h',
            'hR': 'D_3d',
            'hP': 'D_6h',
            'cP': 'O_h', 'cI': 'O_h', 'cF': 'O_h'})
    csym = OrderedDict({
                'aP': '1  triclinic   primitive      -1     C_i \n',
                'mP': '2  monoclinic  primitive      2/m    C_2h\n',
                'mS': '3  monoclinic  primitive      2/m    C_2h\n',
                'oP': '4  orthorombic primitive      mmm    D_2h\n',
                'oS': '5  orthorombic body-centered  mmm    D_2h\n',
                'oI': '6  orthorombic body-centered  mmm    D_2h\n',
                'oF': '7  orthorombic face-centered  mmm    D_2h\n',
                'tP': '8  tetragonal  primitive      4/mmm  D_4h\n',
                'tI': '9  tetragonal  body-centered  4/mmm  D_4h\n',
                'hR': '10 trigonal    primitive      -3m    D_3d\n',
                'hP': '11 hexagonal   primitive      6/mmm  D_6h\n',
                'cP': '12 cubic       primitive      m3m    O_h \n',
                'cI': '13 cubic       face-centered  m3m    O_h \b',
                'cF': '14 cubic       body-centered  m3m    O_h \n'})
    # international tables numbers -> A. Perlovs numbering
    Number2AP = {
                    1 : 1,
                    2 : 2,
                    3 : 3,
                    4 : 5,
                    5 : 7,
                    6 :  13,
                    7 :  15,
                    8 :  21,
                    9 :  27,
                    10 : 39,
                    11 : 41,
                    12 : 43,
                    13 : 50,
                    14 : 55,
                    15 : 61,
                    16 : 73,
                    17 : 74,
                    18 : 77,
                    19 : 80,
                    20 : 81,
                    21 : 84,
                    22 : 87,
                    23 : 88,
                    24 : 89,
                    25 : 90,
                    26 : 93,
                    27 : 99,
                    28 : 102,
                    29 : 108,
                    30 : 114,
                    31 : 120,
                    32 : 126,
                    33 : 129,
                    34 : 135,
                    35 : 138,
                    36 : 145,
                    37 : 147,
                    38 : 150,
                    39 : 156,
                    40 : 162,
                    41 : 168,
                    42 : 174,
                    43 : 177,
                    44 : 180,
                    45 : 183,
                    46 : 186,
                    47 : 192,
                    48 : 193,
                    49 : 195,
                    50 : 198,
                    51 : 204,
                    52 : 210,
                    53 : 216,
                    54 : 222,
                    55 : 228,
                    56 : 231,
                    57 : 234,
                    58 : 240,
                    59 : 247,
                    60 : 249,
                    61 : 255,
                    62 : 257,
                    63 : 263,
                    64 : 269,
                    65 : 275,
                    66 : 278,
                    67 : 281,
                    68 : 287,
                    69 : 299,
                    70 : 300,
                    71 : 302,
                    72 : 303,
                    73 : 306,
                    74 : 308,
                    75 : 314,
                    76 : 315,
                    77 : 316,
                    78 : 317,
                    79 : 318,
                    80 : 319,
                    81 : 320,
                    82 : 321,
                    83 : 322,
                    84 : 323,
                    85 : 324,
                    86 : 326,
                    87 : 328,
                    88 : 329,
                    89 : 331,
                    90 : 332,
                    91 : 333,
                    92 : 334,
                    93 : 335,
                    94 : 336,
                    95 : 337,
                    96 : 338,
                    97 : 339,
                    98 : 340,
                    99 : 341,
                    100 : 342,
                    101 : 343,
                    102 : 344,
                    103 : 345,
                    104 : 346,
                    105 : 347,
                    106 : 348,
                    107 : 349,
                    108 : 350,
                    109 : 351,
                    110 : 352,
                    111 : 353,
                    112 : 354,
                    113 : 355,
                    114 : 356,
                    115 : 357,
                    116 : 358,
                    117 : 359,
                    118 : 360,
                    119 : 361,
                    120 : 362,
                    121 : 363,
                    122 : 364,
                    123 : 365,
                    124 : 366,
                    125 : 368,
                    126 : 370,
                    127 : 371,
                    128 : 372,
                    129 : 374,
                    130 : 376,
                    131 : 377,
                    132 : 378,
                    133 : 380,
                    134 : 382,
                    135 : 383,
                    136 : 384,
                    137 : 386,
                    138 : 388,
                    139 : 389,
                    140 : 390,
                    141 : 392,
                    142 : 394,
                    143 : 396,
                    144 : 398,
                    145 : 400,
                    146 : 401,
                    147 : 404,
                    148 : 405,
                    149 : 407,
                    150 : 408,
                    151 : 409,
                    152 : 410,
                    153 : 411,
                    154 : 412,
                    155 : 414,
                    156 : 415,
                    157 : 416,
                    158 : 417,
                    159 : 418,
                    160 : 420,
                    160 : 419,
                    161 : 422,
                    162 : 423,
                    163 : 424,
                    164 : 425,
                    165 : 426,
                    166 : 428,
                    167 : 430,
                    168 : 431,
                    169 : 432,
                    170 : 433,
                    171 : 434,
                    172 : 435,
                    173 : 436,
                    174 : 437,
                    175 : 438,
                    176 : 439,
                    177 : 440,
                    178 : 441,
                    179 : 442,
                    180 : 443,
                    181 : 444,
                    182 : 445,
                    183 : 446,
                    184 : 447,
                    185 : 448,
                    186 : 449,
                    187 : 450,
                    188 : 451,
                    189 : 452,
                    190 : 453,
                    191 : 454,
                    192 : 455,
                    193 : 456,
                    194 : 457,
                    195 : 458,
                    196 : 459,
                    197 : 460,
                    198 : 461,
                    199 : 462,
                    200 : 463,
                    201 : 465,
                    202 : 466,
                    203 : 468,
                    204 : 469,
                    205 : 470,
                    206 : 471,
                    207 : 472,
                    208 : 473,
                    209 : 474,
                    210 : 475,
                    211 : 476,
                    212 : 477,
                    213 : 478,
                    214 : 479,
                    215 : 480,
                    216 : 481,
                    217 : 482,
                    218 : 483,
                    219 : 484,
                    220 : 485,
                    221 : 486,
                    222 : 488,
                    223 : 489,
                    224 : 491,
                    225 : 492,
                    226 : 493,
                    227 : 494,
                    229 : 498,
                    230 : 499
                    }

    def __init__(self, atoms):
        cell = atoms.get_cell()
        bl = cell.get_bravais_lattice()
        sg = get_spacegroup(atoms)
        ps = bl.pearson_symbol

        self.bravais_number = list(self.sym.keys()).index(ps) + 1
        self.schoenflies_symbol = self.sym[ps]
        self.crys_sys=self.csym[ps]
        self.basis = 0
        self.sgno=sg.no
        self.apno=self.Number2AP[sg.no]


        self.bravais = cell.get_bravais_lattice()
        self.alat = bl.a / Bohr
        self.boa = cell.cellpar()[1] / cell.cellpar()[0]
        self.coa = cell.cellpar()[2] / cell.cellpar()[0]
        self.rbas = sg.scaled_primitive_cell
        self.blat = self.boa*self.alat
        self.clat = self.coa*self.alat
        self.alpha=cell.cellpar()[3]
        self.beta=cell.cellpar()[4]
        self.gamma=cell.cellpar()[5]


    def __repr__(self):
        desc = self.bravais.longname
        s = "<{}({})>".format(self.__class__.__name__, desc)
        return s


class SiteData(object):

    def __init__(self, coords, occupancy=[]):
        self.x, self.y, self.z = coords
        self.occupancy = occupancy
        self.noq = len(occupancy)
        self.irefq=1
        self.imq=1
        self.qmtet = 0.
        self.qmphi = 0.
#        self.__class__.SITES_LIST.append(self)
        self.id = 1 #len(self.__class__.SITES_LIST)
        self.itoq=0

    def __repr__(self):
        s = "<SiteData(id={:d}, noq={:d})>".format(self.id, self.noq)
        return s


class TypeData(object):

    def __init__(self, symbol, concentration=1., iq=None,site=None, mesh=None,iqat=None):
        self.symbol = symbol
        self.atomic_number = Atom(symbol).number
        self.nat = 1
        self.concentration = concentration
        self.im = 1
        self.site = site
        self.iqat = iqat
        self.id = iq

    def __repr__(self):
        s = "<TypeData(\'{}\', id={:d}, concentration={:.4f}, site={:d},iqat={})>".format(
            self.symbol, self.id, self.concentration, self.site.id,self.iqat)
        return s


class ExponentialMeshData(object):
    def __init__(self):
        self.id = 1
        self.type = "exponential"
        self.rws = 0.0
        self.rmt = 0.0
        self.jws = 721
        self.jmt = 0
        self.r_1 = 1e-6
        self.dx  = 2e-2
    def __repr__(self):
        s = "<{}(id={:d}, type=\"{}\")>".format(self.__class__.__name, self.id,
             self.type)
        return s


class PotFile(object):
    VERSION_STRING = "6  (21.05.2007)"

    def __init__(self, atoms, filename=None, title=None, system=None,sysfilename=None,directory=None):
        self.filename = filename
        self.directory = directory
        self.atoms = atoms
        self.title = title
        self.system = system
        if system is None:
            self.system = atoms.get_chemical_formula()
        self.package='SPRKKR'
        self.sysfilename=sysfilename


        self.nm = 1
        self.nq = 1
        self.nt = 2
        self.ns = 2
        self.ne = 30
        self.irel = 3
        self.ibzint = 2
        self.nktab = 0
        self.scf_status='START'
        self.info = None
        self.xc_pot = "VWM"
        self.scf_alg = "BROYDEN2"
        self.fullpot=False
        self.nonmag=False
        self.semicore=False
        self.lloyd='F'
        self.scf_iter = 0
        self.scf_mix = 2e-1
        self.scf_tol = 1e-5
        self.breitint = False
        self.orbpol = None
        self.extfield = False
        self.blcoupl = False
        self.bext = 0
        self.kmrot = 0
        self.rmsavv=999999.
        self.rmsavb=999999.

        self.qmvec = np.array([0,0,0])
        self.ef = 999999.
        self.vmtz=0.7
        self.cartesian=True
        self.basscale = np.array([1.0,1.0,1.0])

        self.vref=4.0000
        self.rmtref=0.00000

        self.tabncore=[0,0,0,2,2,2,2,2,2,2,2,10,10,10,10,10,10,10,10,18,18,
                       18,18,18,18,18,18,18,18,18,28,28,28,28,28,28,28,36,
                       36,36,36,36,36,36,36,36,36,36,46,46,46,46,46,46,46,
                       54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,
                       68,68,68,68,68,68,68,68,78,78,78,78,78,78,78,86,86,
                       86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,86,
                       86,86,86,86,86]

        # create lattice, mesh, site and types from ase.atoms
        self.ld = LatticeData(atoms)
        self.all_at = []
        self.all_sd = []
        self.all_nm = []

        atoms = self.atoms
#Use symmetry in order to find equivalet sites
        sg = get_spacegroup(atoms)
        unique_sites=sg.unique_sites(atoms.get_scaled_positions())
        all_sites = atoms.get_scaled_positions()
        occupancy = get_occupancy(atoms)
        type_table=np.zeros(shape=(len(unique_sites),len(all_sites)),dtype=int)
        #type_table [ unique_sites x all_sites ]: 1 if the site belongs to theequivalence class
        #iqref [ all_sites ] = equivalence class number
        #iqat { unique_site_no } = [ all_sites of given class ]

        tol=1e-5
        iqat={}
        iqref=[]
        itoq={}

        for i in range(len(all_sites)):
            iqref.append(-1)

        if len(unique_sites)!=len(all_sites):
            LOGGER.info(f"==========================================")
            LOGGER.info(f"Symmetry adaped list of types will be used")
            LOGGER.info(f"==========================================")

            for ui in range(len(unique_sites)):
                su,k=sg.equivalent_sites(unique_sites[ui])
                lst=[]
                for i in range(len(su)):
                    for j in range(len(all_sites)):
                        diff=np.abs(su[i]-all_sites[j])<=[tol,tol,tol]
                        if diff.all():
                            type_table[ui][j]=1
                            lst.append(j)
                            iqref[j]=ui
                            if occupancy[j]!=occupancy[ui]:
                                LOGGER.warning(f"\n OCCUPANCY NOT CONFORM WITH THE CRYSTAL SYMMETRY \n")
                                LOGGER.warning(f"Occupation of SITE {j:d} will be taken from SITE {ui:d} \n")
                                occupancy[j]=occupancy[ui]
                iqat[ui]=lst
            LOGGER.info(f"Table of equivalent types \n \"{type_table}\"")
        else:
            for j in range(len(all_sites)):
                type_table[j][j]=1
                iqat[j]=[j]
                iqref[j]=j
#
        dsf={} #dictionary connecting unique site with occupation
        it=0
        for ui in range(len(unique_sites)):
            lsf=[]
            for symbol, concentration in occupancy[ui].items():
                lsf.append(it)
                it=it+1
            dsf[ui]=lsf

        for i, sitexyz in enumerate(all_sites):
            ui=iqref[i]
            sd = SiteData(sitexyz)
            sd.noq=len(occupancy[ui])
            sd.itoq=ui
            sd.irefq=ui
            sd.imq=ui

            ll=0
            lsf=dsf[ui]
            for symbol, concentration in occupancy[ui].items():
                td = TypeData(symbol, concentration,iq=i, site=lsf[ll],iqat=iqat[ui])
                sd.occupancy.append(td)
                self.all_at.append(td)
                ll=ll+1
            self.all_sd.append(sd)
            self.all_sd[-1].id=len(self.all_sd)

        # set the mesh to the atom types
        for ui in range(len(unique_sites)):
            md = ExponentialMeshData()
            md.id=ui
            self.all_nm.append(md)

        key_max = max(dsf.keys(), key=(lambda k: dsf[k]))
        self.nq=len(all_sites)
        self.nt=max(dsf[key_max])+1
        self.nm=len(self.all_nm)

    def write(self):
        pot=os.path.join(self.directory, self.filename)
        with open(pot, "w") as fd:
            fd.write(str(self))

    def write_sys(self):
         sys=os.path.join(self.directory, self.sysfilename)
         with open(sys, "w") as fd:
            fd.write(str(self.sysfile()))

    def sysfile(self):
        filestring = "system data-file created by python ase2sprkkr \n"
        filestring += self.sysfilename+"\n"
        filestring += "xband-version\n"
        filestring += "5.0\n"
        # It would be really cool to support lower dimensions...one day.
        filestring += "dimension\n"
        filestring += "3D\n"
        filestring += "Bravais lattice\n"
        filestring += self.ld.crys_sys
        filestring += "space group number (ITXC and AP)\n"
        filestring += "%5i%5i"%(self.ld.sgno,self.ld.apno)+"\n"
        filestring += "structure type\n"
        filestring += "UNKNOWN\n"
        filestring += "lattice parameter A  [a.u.]\n"
        filestring += "%18.12f\n"%self.ld.alat
        filestring += "ratio of lattice parameters  b/a  c/a\n"
        filestring += "%18.12f%18.12f\n"%(self.ld.boa,self.ld.coa)
        filestring += "lattice parameters  a b c  [a.u.]\n"
        filestring += "%18.12f%18.12f%18.12f\n"%(self.ld.alat,self.ld.blat,self.ld.clat)
        filestring += "lattice angles  alpha beta gamma  [deg]\n"
        filestring += "%18.12f%18.12f%18.12f\n"%(self.ld.alpha,self.ld.beta,self.ld.gamma)
        filestring += "primitive vectors     (cart. coord.) [A]\n"
        for vec in self.ld.rbas:
            for p in vec:
                filestring += "%18.12f"%p
            filestring += "\n"
        # Get number of sites and fill out with empty spheres if the sites are not fully filled
        filestring += "number of sites NQ\n"
        filestring += "%3i\n"%self.nq
        filestring += " IQ ICL     basis vectors     (cart. coord.) [A]                      RWS [a.u.]  NLQ  NOQ ITOQ\n"
        s=""
        rws=0.0
        angmom=4
        for sd in self.all_sd:
            filestring += "%3i%4i%18.12f%18.12f%18.12f  %18.12f%4i%5i "%(sd.id,sd.irefq+1,sd.x,sd.y,sd.z,rws,angmom,sd.noq)
            for at in sd.occupancy:
                if at.id+1 == sd.id:
                   filestring += "%3i"%(at.site+1)
            filestring+="\n"
        filestring+="number of sites classes NCL \n"
        filestring += "%3i\n"%(self.nt)
        filestring+="ICL WYCK NQCL IQECL (equivalent sites)\n"
        itdone= [False for i in range(self.nt)]
        for iat, at in enumerate(self.all_at):
            if (not itdone[at.site]):
                filestring += "%3i   %1s%5i"%(at.site+1,'-',len(at.iqat))
                for key in at.iqat:
                    filestring += "%3i"%(key+1)
                filestring+="\n"
                itdone[at.site]=True

        filestring += "number of atom types NT\n"
        filestring += "%3i\n"%self.nt

        filestring += " IT  ZT  TXTT  NAT  CONC  IQAT (sites occupied)\n"
        itdone= [False for i in range(self.nt)]
        for iat, at in enumerate(self.all_at):
            if (not itdone[at.site]):
                filestring += " %2i%4i  %8s%5i%6.3f"%(at.site+1,at.atomic_number,at.symbol,len(at.iqat),at.concentration)
                for key in at.iqat:
                    filestring += "%3i"%(key+1)
                filestring+="\n"
                itdone[at.site]=True
        # Average Wigner-Seitz radi
# =============================================================================
#         rws = 0.
#         iq = 0
#         icl = 0
#         itoq = 0
#         for a in self.cell.atomdata:
#             icl += 1
#             itoqs = []
#             for sp in a[0].species:
#                 itoq += 1
#                 itoqs.append(itoq)
#             for b in a:
#                 iq += 1
#                 if self.minangmom:
#                     angmom = max(max([ed.angularmomentum[ed.elementblock[spcs]] for spcs in b.species])+1,self.minangmom)
#                 else:
#                     angmom = max([ed.angularmomentum[ed.elementblock[spcs]] for spcs in b.species])+1
#                 v = mvmult3(self.cell.latticevectors,b.position)
#                 filestring += "%3i%4i%18.12f%18.12f%18.12f  %18.12f%4i%5i "%(iq,icl,v[0],v[1],v[2],rws,angmom,len(a[0].species))
#                 for i in itoqs:
#                     filestring += "%3i"%i
#                 filestring += "\n"
#         filestring += "number of sites classes NCL\n"
#         filestring += "%3i\n"%len(self.nt)
#         filestring += "ICL WYCK NQCL IQECL (equivalent sites)\n"
#         iq = 0
#         icl = 0
#         for a in self.cell.atomdata:
#             icl += 1
#             filestring += "%3i   %1s%5i"%(icl,'-',len(a))
#             for b in a:
#                 iq += 1
#                 filestring += "%3i"%iq
#             filestring += "\n"
#         filestring += "number of atom types NT\n"
#         nt = 0
#         for a in self.cell.atomdata:
#             nt += len(a[0].species)
#         filestring += "%3i\n"%nt
#         filestring += " IT  ZT  TXTT  NAT  CONC  IQAT (sites occupied)\n"
#         iq = 0
#         it = 0
#         for a in self.cell.atomdata:
#             corr = 0
#             for sp,conc in a[0].species.iteritems():
#                 it += 1
#                 filestring += " %2i%4i  %8s%5i%6.3f"%(it,ed.elementnr[sp],sp,len(a),conc)
#                 iq -= corr*len(a)
#                 for b in a:
#                     iq += 1
#                     filestring += "%3i"%iq
#                 corr = 1
#                 filestring += "\n"
# =============================================================================
        return filestring



    def __str__(self):
        separator = "*"*80 + "\n"
        #**************************************************
        s  = separator
        s += 'HEADER    SCF-start data created by ase2sprkkr'
        s +=  datetime.now().ctime() + "\n"
        #**************************************************
        s += separator
        #**************************************************
        s += "{:<10s}{}\n".format("TITLE", self.title)
        s += "{:<10s}{}\n".format("SYSTEM", self.system)
        s += "{:<10s}{}\n".format("PACKAGE", self.package)
        s += "{:<10}{}\n".format("FORMAT", self.VERSION_STRING)
        #**************************************************
        s += separator
        #**************************************************
        s += "{}\n".format("GLOBAL SYSTEM PARAMETER")
        s += "{:<18}{:>2d}\n".format("NQ", self.nq)
        s += "{:<18}{:>2d}\n".format("NT", self.nt)
        s += "{:<18}{:>2d}\n".format("NM", self.nm)
        s += "{:<18}{:>2d}\n".format("IREL", self.irel)
        #**************************************************
        s += separator
        #**************************************************
        s += "{}\n".format("SCF-INFO")
        s += "{:<10}{}\n".format("INFO", "NONE" if self.info is None else self.info)
        s += "{:<10}{}\n".format("SCFSTATUS",self.scf_status)
        s += "{:<10}{}\n".format("FULLPOT", "F" if not self.fullpot else "T")
        s += "{:<12}{}\n".format("BREITINT", "F" if not self.breitint else "T")
        s += "{:<10}{}\n".format("NONMAG","F" if not self.nonmag else "T")
        s += "{:<10}{}\n".format("ORBPOL", "NONE" if self.orbpol is None else
                                self.orbpol)
        s += "{:<12}{}\n".format("EXTFIELD", "F" if not self.extfield else "T")
        s += "{:<12}{}\n".format("BLCOUPL", "F" if not self.blcoupl else "T")
        s += "{:<18}{:1.10f}\n".format("BEXT", self.bext)
        s += "{:<10}{}\n".format("SEMICORE",self.semicore)
        s += "{:<10}{}\n".format("LLOYD",self.lloyd)
        s += "{:<18}{:>2d}\n".format("NE", self.ne)
        s += "{:<18}{:>2d}\n".format("IBZINT", self.ibzint)
        s += "{:<18}{:>2d}\n".format("NKTAB", self.nktab)
        s += "{:<10}{}\n".format("XC-POT", self.xc_pot)
        s += "{:<10}{}\n".format("SCF-ALG", self.scf_alg)
        s += "{:<18}{:d}\n".format("SCF-ITER", self.scf_iter)
        s += "{:<18}{:1.10f}\n".format("SCF-MIX", self.scf_mix)
        s += "{:<18}{:1.10f}\n".format("SCF-TOL", self.scf_tol)
        s += "{:<18}{:1.10f}\n".format("RMSAVV", self.rmsavv)
        s += "{:<18}{:1.10f}\n".format("RMSAVB", self.rmsavb)
        s += "{:<10}{:>20.10f}\n".format("EF", self.ef)
        s += "{:<10}{:1.10f}\n".format("VMTZ", self.vmtz)
        s += separator
        #**************************************************
        s += "{}\n".format("LATTICE")
        s += "{}\n".format("SYSDIM    3D")
        s += "{}\n".format("SYSTYPE   BULK")
        desc = self.ld.bravais.longname
        pg = ""
        s += "{:<18}{:>2d}  {:>31}  {}  {}\n".format("BRAVAIS",
             self.ld.bravais_number, desc, pg, self.ld.schoenflies_symbol)
        s += "{:<12}{:1.10f}\n".format("ALAT", self.ld.alat)
        for i, row in enumerate(self.ld.rbas):
            s += ("A({:d})" + "{:>20.10f}" * 3 + "\n").format(i+1, *row)
        s += separator
        #**************************************************
        s += "{}\n".format("SITES")
        s += "{:<12}{}\n".format("CARTESIAN", "F" if not self.cartesian else "T")
        s += ("{:<10}" + "{:>15.10f}"*3 + "\n").format("BASSCALE", *self.basscale)
        s += "{}\n".format("        IQ      QX              QY              QZ")
        for sd in self.all_sd:
            s += ("{:>10d}"+"{:>15.10f}"*3+"\n").format(sd.id, sd.x, sd.y, sd.z)
        s += separator
        #**************************************************
        s += "{}\n".format("OCCUPATION")
        s += "{}\n".format("        IQ     IREFQ       IMQ       NOQ  ITOQ  CONC")
        for sd in self.all_sd:
            s += ("{:>10d}"+
                  "{:>10d}"*3).format(sd.id, sd.irefq+1,sd.imq+1,sd.noq)
            for at in sd.occupancy:
                if at.id+1 == sd.id:
                    s += "{:>3d} {:1.3f}".format(at.site+1,at.concentration)
            s += "\n"
        s += separator
        #**************************************************
        s += "{}\n".format("REFERENCE SYSTEM")
        s += "{:<18}{:>2d}\n".format("NREF", self.nm)
        s += "{}\n".format("      IREF      VREF            RMTREF")
        for iref in range(self.nm):
            s += ("{:>10d}"+"{:>15.10f}"*2+"\n").format(iref+1, self.vref, self.rmtref)
        s += separator
        #**************************************************
        s += "{}\n".format("MAGNETISATION DIRECTION")
        s += "{:<18}{:>2d}\n".format("KMROT", self.kmrot)
        s += ("{:<10}" + "{:>15.10f}"*3 + "\n").format("QMVEC", *self.qmvec)
        s += "{}\n".format("         IQ      QMTET           QMPHI")
        for sd in self.all_sd:
            s += "{:>10d} {:>15.10f} {:15.10f}\n".format(sd.id,
                  sd.qmtet, sd.qmphi)
        s += separator
        #**************************************************
        s += "{}\n".format("MESH INFORMATION")
        s += "MESH-TYPE {}\n".format(self.all_nm[0].type.upper())
        s += "{}\n".format("   IM      R(1)            DX         JRMT      RMT        JRWS      RWS")
        for md in self.all_nm:
            s += "{:>5d}".format(md.id+1)
            s += "{:>18.10f}".format(md.r_1)
            s += "{:>18.10f}".format(md.dx)
            s += "{:>5d}".format(md.jmt)
            s += "{:>18.10f}".format(md.rmt)
            s += "{:>5d}".format(md.jws)
            s += "{:>18.10f}\n".format(md.rws)
        s += separator
        #**************************************************
        s += "{}\n".format("TYPES")
        s += "{}\n".format("   IT     TXTT        ZT     NCORT     NVALT    NSEMCORSHLT")

        itdone= [False for i in range(self.nt)]
        for iat, at in enumerate(self.all_at):
            if (not itdone[at.site]):
                s += "{:>5d}".format(at.site+1)
                s += "{:>9}".format(at.symbol)
                s += "{:>10d}".format(at.atomic_number)
                s += "{:>10d}".format(self.tabncore[at.atomic_number])
                s += "{:>10d}".format(at.atomic_number-self.tabncore[at.atomic_number])
                s += "{:>10d}\n".format(0)
                itdone[at.site]=True

        s += separator

        return s
