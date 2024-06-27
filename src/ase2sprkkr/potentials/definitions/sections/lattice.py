import numpy as np

from ase.units import Bohr
from ase.cell import Cell

from ...potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from ...potential_sections import PotentialSection

from ....ase.pbc import check_symmetry
from ....common.grammar_types import DefKeyword, Sequence, Table, Integer, Array
from ....physics.lattice_data import Pearson
from ....sprkkr.atoms_region import AtomsRegion


class LatticeSection(PotentialSection):

  """ This section retrieves the lattice geometry and
      it creates (during reading) the ASE Cell object """

  def _set_from_atoms(self, atoms, write_io_data):
      pbc = atoms.pbc.sum()
      if pbc == 3:
          self.SYSDIM = '3D'
          bcell = atoms.cell
      elif pbc == 2:
          if not ('left' in atoms.regions and 'right' in atoms.regions and 'central' in atoms.regions):
              raise ValueError("To run a 2D calculation, an atoms object have to have defined "
              "'left', 'right' and 'central' region.")
          for i in ('left', 'right'):
              check_symmetry(atoms.regions[i].pbc, True, lambda x: f'{i.capitalize()} region of a 2D calculations have to be symmetric in all dimensions')
          check_symmetry(atoms.regions['central'].pbc, [True, True, False], lambda x: 'Central region of a 2D calculations have to be symmetric just in X and Y axis.')
          check_symmetry(atoms.pbc, [True, True, False], lambda x: 'For a 2D calculation, atoms have to be periodic in X and Y axes.')
          self.SYSDIM = '2D'
          bcell = atoms.regions['left'].cell
      else:
          raise ValueError(f"I don't know which calculation I should run with the periodicity of type: {atoms.pbc}.")

      bravais_lattice = bcell.get_bravais_lattice()
      alat = bravais_lattice.a
      pearson_symbol = bravais_lattice.pearson_symbol
      self['BRAVAIS'].set(Pearson.from_symbol(pearson_symbol).xband_data())
      if self.SYSDIM() == '3D':
         self.SYSTYPE = 'BULK'
         write_io_data['sites_order'] = slice(None)
         self.A_L3.clear()
         self.A_R3.clear()
         self.NQ_L.clear()
         self.NQ_R.clear()
         cell = atoms.cell
      elif self.SYSDIM() == '2D':
         for i in 'left', 'right', 'central':
             if not i in atoms.regions:
                  raise ValueError("For a 2D problem, the 'left', 'right' and 'central' regions have to be defined")
         self.SYSTYPE = 'LIV' if atoms.regions['right'].only_vacuum_atoms() else 'LIR'
         cell = atoms.cell
         if (cell[2] == [0,0,0]).all():
             if (atoms.regions['central'] == [0,0,0]).all():
                  raise ValueError("For a 2D problem, the third (the z-) cell vector has to be specified either for "
                                   " atoms object, or for its central region")
             cell[2] = atoms.regions['left'].cell[2] + atoms.regions['right'].cell[2] + atoms.regions['central'].cell[2]
         self.A_L3 = atoms.regions['left'].cell[2] / alat
         self.A_R3 = atoms.regions['right'].cell[2] / alat
         self.NQ_L = len(atoms.regions['left'])
         self.NQ_R = len(atoms.regions['right'])
         for i,j in [('left', 'right'), ('left', 'central'), ('right', 'central')]:
             if atoms.regions[i].shared_ids_with(atoms.regions[j]):
                 raise ValueError(f"{i.capitalize()} and {j} region can not share the same site")
         la = len(atoms)
         if self.NQ_L() + self.NQ_R() + len(atoms.regions['central']) != la:
             raise ValueError("The left, central and right do not cover all the atoms in the sample")

         def order(region):
             region = atoms.regions[region]
             z = region.positions[:,2]
             return region.ids[np.argsort(z)]

         write_io_data['sites_order'] = np.concatenate(list(map(order,['left', 'central', 'right'])))
         for i,j in zip(range(len(atoms)), write_io_data['sites_order']):
             if i!=j:
                break
         else:
             # optimalization - do not reorder
             write_io_data['sites_order'] = slice(None)
      else:
         raise ValueError(f'{self.SISDIM()} problem periodicity type not implemented')
      write_io_data['lattice.alat'] = alat
      self['ALAT'].set(alat / Bohr)
      self['SCALED_PRIMITIVE_CELL'].set(cell / alat)

  def _update_atoms(self, atoms, read_io_data):
      alat = self['ALAT']() * Bohr
      cell = Cell(self['SCALED_PRIMITIVE_CELL']() * alat)
      read_io_data['lattice.alat'] = alat
      if self.SYSDIM() == '3D':
          regions = []
          pbc = [True, True, True]
      elif self.SYSDIM() == '2D':
          lc = cell.copy()
          lc[2] = self.A_L3() * alat
          rc = cell.copy()
          rc[2] = self.A_R3() * alat
          cc = cell.copy()
          cc[2] -= lc[2] + rc[2]
          regions = [
              AtomsRegion('left',   slice(None,self.NQ_L()),            lc, [True, True, True]),       # NOQA E241
              AtomsRegion('central',slice(self.NQ_L(), -self.NQ_R()),   cc, [True, True, False]),      # NOQA E241
              AtomsRegion('right',  slice(-self.NQ_R(), None),          rc, [True, True, True])        # NOQA E241
          ]
          pbc = [True, True, False]
          cell[2]+=lc[2] + rc[2]
      else:
          raise ValueError(f'Unknown problem periodicity type {self.SYSDIM()}')

      def update(atoms):
          atoms.cell = cell
          atoms.pbc = pbc
          atoms.set_regions(regions)

      read_io_data.apply_on_atoms(update, atoms)


class LatticeSectionDefinition(PotSectionDefinition):

  def __init__(self, name='LATTICE', **kwargs):
      V = PotValueDefinition
      members = [
          V('SYSDIM', DefKeyword('3D', '2D')),
          V('SYSTYPE', DefKeyword('BULK', 'LIV', 'LIR')),
          V('BRAVAIS', Sequence(int, str, str, str, str, allowed_values = (i.xband_data() for i in Pearson.pearsons.values()))),
          V('ALAT', float),
          # Keywords and thus the numbering has just (or at least) 10 char long
          V('SCALED_PRIMITIVE_CELL', Table([float] * 3, numbering=Integer(prefix='A(', postfix=')', after_format='<10'),length=3, format='>22.14f')),
          V('NQ_L', int, is_optional=True),
          V('A_L(3)', Array(float, length=3, format='<10'), is_optional=True),
          V('NQ_R', int, is_optional=True),
          V('A_R(3)', Array(float, length=3, format='<10'), is_optional=True),
      ]
      super().__init__(name, members, has_hidden_members=True)

  def validate(self, data, why:str='set'):
      d3 = data['SYSDIM'] == '3D'
      if d3 != (data['SYSTYPE'] == 'BULK'):
          raise ValueError('LATTICE.SYSDIM have to be 2D to create LIV or LIR structures')
      op = (lambda x: x is not None) if d3 else (lambda x: x is None)
      if op(data['NQ_L']):
          raise ValueError('LATTICE.NQ_L can be specified only for 2D structures')
      if op(data['A_L(3)']):
          raise ValueError('LATTICE.A_L(3) can be specified only for 2D structures')
      if op(data['NQ_R']):
          raise ValueError('LATTICE.NQ_R can be specified only for 2D structures')
      if op(data['A_R(3)']):
          raise ValueError('LATTICE.A_R(3) can be specified only for 2D structures')
      return True

  result_class = LatticeSection


section = LatticeSectionDefinition
