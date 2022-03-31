"""
This module contains LatticeData object, that provides various informations about a lattice
"""

from ase.spacegroup import get_spacegroup
from ase.units import Bohr

class LatticeData(object):
    """
    This class provide various informations about the lattice
    """

    #: Mapping between two symmetries notation
    symmetries = {
            'aP': 'C_i',
            'mP': 'C_2h', 'mS': 'C_2h',
            'oP': 'D_2h', 'oS': 'D_2h', 'oI': 'D_2h', 'oF': 'D_2h',
            'tP': 'D_4h', 'tI': 'D_4h',
            'hR': 'D_3d',
            'hP': 'D_6h',
            'cP': 'O_h', 'cI': 'O_h', 'cF': 'O_h'
    }

    #: Mapping between Pearson symbols and properties of the cell
    cell_symmetries = {
        'aP': (1,  'triclinic',   'primitive',      '-1',     'C_i'),
        'mP': (2,  'monoclinic',  'primitive',      '2/m',    'C_2h'),
        'mS': (3,  'monoclinic',  'primitive',      '2/m',    'C_2h'),
        'oP': (4,  'orthorombic', 'primitive',      'mmm',    'D_2h'),
        'oS': (5,  'orthorombic', 'body-centered',  'mmm',    'D_2h'),
        'oI': (6,  'orthorombic', 'body-centered',  'mmm',    'D_2h'),
        'oF': (7,  'orthorombic', 'face-centered',  'mmm',    'D_2h'),
        'tP': (8,  'tetragonal',  'primitive',      '4/mmm',  'D_4h'),
        'tI': (9,  'tetragonal',  'body-centered',  '4/mmm',  'D_4h'),
        'hR': (10, 'trigonal',    'primitive',      '-3m',    'D_3d'),
        'hP': (11, 'hexagonal',   'primitive',      '6/mmm',  'D_6h'),
        'cP': (12, 'cubic',       'primitive',      'm3m',    'O_h'),
        'cF': (13, 'cubic',       'face-centered',  'm3m',    'O_h'),
        'cI': (14, 'cubic',       'body-centered',  'm3m',    'O_h')
    }

    #: Translation of international tables numbers to A. Perlovs numbering.
    international_numbers_to_AP = {
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

    @staticmethod
    def bravais_number_from_pearson_symbol(pearson_symbol):
        """ Return lattice data for the given pearson symbol 
        
        Returns
        -------
        lattice_data: tuple
            See `:class:LatticeData.cell_symmetries<ase2sprkkr.physics.lattice_data>`
        """
        return list(LatticeData.symmetries.keys()).index(pearson_symbol) + 1

    def __init__(self, atoms):
        cell = atoms.get_cell()
        bl = cell.get_bravais_lattice()
        sg = get_spacegroup(atoms)
        ps = bl.pearson_symbol

        self.bravais_number = self.bravais_number_from_pearson_symbol(ps)
        self.schoenflies_symbol = self.symmetries[ps]
        self.crys_sys=self.cell_symmetries[ps]
        self.basis = 0
        self.sgno=sg.no
        self.apno=self.international_numbers_to_AP[sg.no]

        self.bravais = cell.get_bravais_lattice()

        self.boa = cell.cellpar()[1] / cell.cellpar()[0]
        self.coa = cell.cellpar()[2] / cell.cellpar()[0]
        self.rbas = sg.scaled_primitive_cell

        self.alat = bl.a / Bohr
        self.blat = self.boa*self.alat
        self.clat = self.coa*self.alat

        self.alpha=cell.cellpar()[3]
        self.beta=cell.cellpar()[4]
        self.gamma=cell.cellpar()[5]
