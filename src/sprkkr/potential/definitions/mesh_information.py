from ..potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from ..potential_sections import UniqueListSection

from  ...common.grammar_types import DefKeyword, Array, Table, Integer
from ...sprkkr.radial_meshes import ExponentialMesh

class MeshInformationSection(UniqueListSection):
  _value_name = 'meshes'
  _value_class = ExponentialMesh

class MeshInformationSectionDefinition(PotSectionDefinition):

  def __init__(self, name='MESH INFORMATION', **kwargs):
      V = PotValueDefinition
      members = [
          V('MESH-TYPE', DefKeyword('EXPONENTIAL')),
          V('DATA', Table({'R(1)': 0., 'DX' :0., 'JRMT': 0, 'RMT': 0., 'JWRS': 0, 'RWS': 0.}, numbering='IM')),
      ]
      super().__init__(name, members, has_hidden_members=True)

  result_class = MeshInformationSection
