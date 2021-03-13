from .potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from .atom_sections import AtomSection

from ase.units import Bohr
from ase.spacegroup import get_spacegroup
from ase.lattice import bravais_classes
from  ..common.grammar_types import DefKeyword, Sequence, Table, Integer

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
    'cI': (13, 'cubic',       'face-centered',  'm3m',    'O_h'),
    'cF': (14, 'cubic',       'body-centered',  'm3m',    'O_h')
}

cell_symmetries_lookup = {v : k for k,v in cell_symmetries.items() }

class LatticeSection(AtomSection):

  _lattice_params = {
        'a' : 'ALAT',
        'b': 'BLAT',
        'c': 'CLAT',
        'alpha' : 'ALPHALAT'
  }

  def _set_from_atoms(self):
      aw = self._atoms_wrapper
      bravais_lattice = aw.bravais_lattice
      pearson_symbol = bravais_lattice.pearson_symbol
      lattice_type = cell_symmetries[pearson_symbol]

      self['BRAVAIS'].set(lattice_type)
      for n,k in self._lattice_params.items():
          val = getattr(bravais_lattice, n, None)
          if val is None:
              self[k].clear()
          else:
              self[k].set(val / Bohr)
      self['SCALED_PRIMITIVE_CELL'].set(aw.scaled_primitive_cell)

  def _process(self):
      try:
        cell = cell_symmetries_lookup[self['BRAVAIS']()]
        bravais_class = bravais_classes[cell]
        args = {}
        for n,k in self._lattice_params.items():
            if self[k]():
               args[n] = self[k]() * Bohr
        cell = bravais_class(**args)
      except KeyError:
        cell = None
      self._atoms_io_data.cell = cell

class LatticeSectionDefinition(PotSectionDefinition):

  def __init__(self, name, **kwargs):
      V = PotValueDefinition
      members = [
          V('SYSDIM', DefKeyword('3D')),
          V('SYSTYPE', DefKeyword('BULK')),
          V('BRAVAIS', Sequence(int, str, str, str, str, allowed_values = cell_symmetries.values())),
          V('ALAT', float)
        ] + [
          V(i, float, is_optional=i!='ALAT') for
               i in LatticeSection._lattice_params.values()
        ] + [
          V('SCALED_PRIMITIVE_CELL', Table([float]*3, numbering=Integer(prefix='A(', postfix=')'),length=3 )),
      ]
      super().__init__(name, members, has_hidden_members=True)

  result_class = LatticeSection
