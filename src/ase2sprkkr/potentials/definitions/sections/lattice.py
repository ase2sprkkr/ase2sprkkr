from ...potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from ...potential_sections import PotentialSection

from ase.units import Bohr
from ase.lattice import bravais_classes
from ase.cell import Cell
from ....common.grammar_types import DefKeyword, Sequence, Table, Integer
from ....physics.lattice_data import LatticeData

class LatticeSection(PotentialSection):

  """ This section retrieves the lattice geometry and
      it creates (during reading) the ASE Cell object """

  def _set_from_atoms(self, atoms, write_io_data):
      bravais_lattice = atoms.cell.get_bravais_lattice()
      pearson_symbol = bravais_lattice.pearson_symbol
      self['BRAVAIS'].set(LatticeData.cell_symmetries[pearson_symbol])
      alat = bravais_lattice.a
      self['ALAT'].set(alat / Bohr)
      self['SCALED_PRIMITIVE_CELL'].set(atoms.cell / alat)
      write_io_data['lattice.alat'] = alat

  def _update_atoms(self, atoms, read_io_data):
      alat = self['ALAT']() * Bohr
      cell = Cell(self['SCALED_PRIMITIVE_CELL']() * alat)
      read_io_data['lattice.cell'] = cell
      read_io_data['lattice.alat'] = alat
      if atoms:
         atoms.cell = cell

class LatticeSectionDefinition(PotSectionDefinition):

  def __init__(self, name='LATTICE', **kwargs):
      V = PotValueDefinition
      members = [
          V('SYSDIM', DefKeyword('3D')),
          V('SYSTYPE', DefKeyword('BULK')),
          V('BRAVAIS', Sequence(int, str, str, str, str, allowed_values = LatticeData.cell_symmetries.values())),
          V('ALAT', float),
          #Keywords and thus the numbering has just (or at least) 10 char long
          V('SCALED_PRIMITIVE_CELL', Table([float]*3, numbering=Integer(prefix='A(', postfix=')', format='<10'),length=3)),
      ]
      super().__init__(name, members, has_hidden_members=True)

  result_class = LatticeSection

section = LatticeSectionDefinition
