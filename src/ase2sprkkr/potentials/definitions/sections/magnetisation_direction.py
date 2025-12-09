from ...potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from ...potential_sections import PotentialSection

from ....common.grammar_types import Table, Array


class MagnetisationDirectionSection(PotentialSection):
    def _depends_on(self):
        return ['TYPES']

    def _set_from_atoms(self, atoms, write_io_data):
        self.DATA_IQ = [ site.magnetisation for site in
                            atoms.sites[write_io_data['sites_order']] ]
        self.DATA_IT = [ type.magnetisation for type in
                write_io_data.types.value_to_class_id ]
        self.KMROT = atoms.info.get('SPRKKR_KMROT', 0)
        self.QMVEC = atoms.info.get('SPRKKR_QMVEC', [0.0,0.0,0.0])

    def _update_atoms(self, atoms, read_io_data):
        atoms.info['SPRKKR_KMROT'] = self.KMROT()
        atoms.info['SPRKKR_QMVEC'] = self.QMVEC()

        for site, d  in zip(atoms.sites, self.DATA_IQ()):
            site.magnetisation = d

        if self.DATA_IT() is not None:
            for typ,d  in zip(read_io_data['types'], self.DATA_IT()):
                typ.magnetisation = d


class MagnetisationDirectionSectionDefinition(PotSectionDefinition):

  array_name = 'charge_moments'

  def __init__(self, name='CHARGE MOMENTS', **kwargs):
      V = PotValueDefinition
      members = [
          V('KMROT', 0),
          V('QMVEC', Array(float, length=3), [0.,0.,0.]),
          V('DATA_IQ', Table({'MTET_Q': float, 'MPHI_Q': float}, numbering='IQ', numbering_format='{:>10}', free_header=True)),
          V('DATA_IT', Table({'MTET_T':float, 'MPHI_T':float, 'MGAM_T': float}, numbering='IT', free_numbering=True, numbering_format='{:>10}',
                                                                                free_header= lambda x: '*' not in x), is_optional=True),
      ]
      super().__init__(name, members, has_hidden_members=True)

  result_class = MagnetisationDirectionSection


section = MagnetisationDirectionSectionDefinition
