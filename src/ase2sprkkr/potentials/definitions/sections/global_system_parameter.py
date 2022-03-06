from ...potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from ...potential_sections import PotentialSection

class GlobalSystemParameter(PotentialSection):

  def _set_from_atoms(self, atoms, write_io_data):
      self['NQ'].set(len(atoms.sites))
      self['NT'].set(write_io_data.types.len_of_unique())
      self['NM'].set(write_io_data.meshes.len_of_unique())

class GlobalSystemParameterDefinition(PotSectionDefinition):

  def __init__(self, name='GLOBAL SYSTEM PARAMETER', **kwargs):
      V = PotValueDefinition
      members = [
        V('NQ', int),
        V('NT', int),
        V('NM', int),
        V('IREL', 3),
        V('NSPIN', int, required = False, is_optional = True)
      ]
      super().__init__(name, members, has_hidden_members=True)

  result_class = GlobalSystemParameter

section = GlobalSystemParameterDefinition
