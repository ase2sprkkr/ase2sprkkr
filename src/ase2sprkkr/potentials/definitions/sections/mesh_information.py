from ...potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from ...potential_sections import UniqueListSection

from ....common.grammar_types import DefKeyword, Table, Array
from ....sprkkr.radial_meshes import ExponentialMesh

from ....sprkkr.io_data import ReadIoData, WriteIoData
from ....sprkkr.sprkkr_atoms import SPRKKRAtoms


class MeshInformationSection(UniqueListSection):
  _value_name = 'meshes'
  _value_class = ExponentialMesh

  def _set_from_atoms(self, atoms:SPRKKRAtoms, write_io_data: WriteIoData):
      ul = getattr(write_io_data, self._value_name)
      self['DATA'].set([ i.to_tuple() for i in ul.iter_unique()])
      fullpot = [ i.fullpot_tuple() for i in ul.iter_unique() ]

      found = None
      for i in fullpot:
          if i is None:
              if found is True:
                  break
              found = False
          else:
              if found is False:
                  break
              found = True
      else:
          if found:
              self['FULLPOT'].set(fullpot)
          return
      raise Exception("FULLPOT mesh information missing for some meshes")

  def _update_atoms(self, atoms:SPRKKRAtoms, read_io_data: ReadIoData):
      super()._update_atoms(atoms, read_io_data)
      if self.FULLPOT() is not None:
          for mesh, full in zip(read_io_data[self._value_name], self['FULLPOT']()):
              mesh.set_fullpot(full)


class MeshInformationSectionDefinition(PotSectionDefinition):

  def __init__(self, name='MESH INFORMATION', **kwargs):
      V = PotValueDefinition
      members = [
          V('MESH-TYPE', DefKeyword('EXPONENTIAL')),
          V('DATA', Table({'R(1)': 1e-6, 'DX' :2e-2, 'JRMT': 0, 'RMT': 0., 'JRWS': 721, 'RWS': 0.},
                    numbering='IM', numbering_format='>5', format='>16',  # header_length=80
                          )
            ),
          V('FULLPOT', Table({'JRNS1': int, 'JRCRI': int, 'NPAN': int, 'JRCUT': Array(int) }, numbering='IM', numbering_format='>5', format='>16'),
            is_optional=True)
      ]
      super().__init__(name, members, has_hidden_members=True)

  result_class = MeshInformationSection


section = MeshInformationSectionDefinition
