"""
Scf info section contains configuration options used in computation
"""

from ...potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from ...potential_sections import PotentialSection

from ....common.grammar_types import DefKeyword, Array, line_string


class ScfInfoSection(PotentialSection):
  """ This section retrieves the atomic positions and
      it creates (during reading) the ASE Atoms object """

  def _depends_on(self):
      return ['SITES']

  def _update_atoms(self, atoms, read_io_data):
      atoms.info['fermi_energy'] = self.EF()
      atoms.info['sprkkr_rmsavv'] = self.RMSAVV()
      atoms.info['sprkkr_rmsavb'] = self.RMSAVB()
      atoms.info['sprkkr_vmtz'] = self.VMTZ()

  def _set_from_atoms(self, atoms, write_io_data):
      info = atoms.info
      if 'fermi_energy' in info:
          self.EF = atoms.info['fermi_energy']
      if 'sprkkr_rmsavv' in info:
          self.RMSAVV = atoms.info['sprkkr_rmsavv']
      if 'sprkkr_rmsavb' in info:
          self.RMSAVB = atoms.info['sprkkr_rmsavb']
      if 'sprkkr_vmtz' in info:
          self.VMTZ = atoms.info['sprkkr_vmtz']
      if self.SCFSTATUS() == 'START':
          for i in atoms.sites:
              if i.potential is None:
                  break
          else:
              self.SCFSTATUS = 'CONVERGED'


class ScfInfoSectionDefinition(PotSectionDefinition):

  def __init__(self, name='SCF-INFO', **kwargs):
      V = PotValueDefinition
      members = [
        V('INFO', line_string, 'NONE'),
        V('SCFSTATUS', DefKeyword('START', 'CONVERGED', 'ITR-BULK')),
        V('FULLPOT', False),
        V('BREITINT', False),
        V('NONMAG', False, alternative_names='NOMAG'),
        V('ORBPOL', str, 'NONE'),
        V('EXTFIELD', False),
        V('BLCOUPL', False),
        V('BEXT', 0.0),
        V('SEMICORE', False),
        V('LLOYD', False),
        V('NE', Array(int), is_optional = True),
        V('IBZINT', int, is_optional = True),
        V('NKTAB', int, is_optional = True),
        V('XC-POT', str, is_optional = True),
        V('SCF-ALG', str, is_optional = True),
        V('SCF-ITER', 0),
        V('SCF-MIX', 0.2),
        V('SCF-TOL', 1e-5),
        V('RMSAVV', 999999.),
        V('RMSAVB', 999999.),
        V('EF', 999999.),
        V('VMTZ', 0.7),
      ]
      super().__init__(name, members, has_hidden_members=True)

  result_class = ScfInfoSection


section = ScfInfoSectionDefinition
