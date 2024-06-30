from ...potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from ...potential_sections import PotentialSection

from ....common.grammar_types import Table, unsigned, Array, Sequence
import numpy as np
from ....sprkkr.sites import Site, Occupation
from ....bindings.spglib import SpacegroupInfo


class OccupationSection(PotentialSection):
  """ This section retrieves the atomic positions and
      it creates (during reading) the ASE Atoms object """

  def _set_from_atoms(self, atoms, write_io_data):
      meshes = write_io_data.meshes
      ref_systems = write_io_data.reference_systems
      types = write_io_data.types

      def onesite(site):
          occ = site.occupation
          gen_occ = ((types.value_to_class_id[i], occ[i]) for i in occ)
          return (
              ref_systems.value_to_class_id[site.reference_system],
              meshes.value_to_class_id[site.mesh],
              len(site.occupation),
              list(gen_occ)
          )

      self['DATA'].set([
        onesite(i) for i in atoms.sites[write_io_data['sites_order']]
      ])

  def _update_atoms(self, atoms, read_io_data):
      data = self['DATA']
      indexes = [ tuple(d) for d in data() ]
      sg_info = SpacegroupInfo.from_atoms(atoms, indexes)
      unique = sg_info.equivalent_sites.mapping

      tags = {}

      def site(i, d):
          ind = unique[i]
          if not ind in tags:
             occ = d['ITOQ CONC']
             occ = Occupation({(read_io_data['types'][i - 1], oc) for i,oc in occ})
             site = Site.create(atoms = atoms,
                         occupation = occ,
                         reference_system = read_io_data['reference_systems'][d['IREFQ'] - 1],
                         mesh = read_io_data['meshes'][d['IMQ'] - 1])
             tags[ind] = site
             return site
          else:
             return Site(tags[ind].site_type)

      atoms.set_sites( np.array([ site(i,d) for i,d in enumerate(data()) ]), sg_info )


class OccupationSectionDefinition(PotSectionDefinition):

  def __init__(self, name='OCCUPATION', **kwargs):
      V = PotValueDefinition

      def validate_row(result):
          return result[2] == len(result[3]) or 'Number of occupations differs from NOQ'

      members = [
          V('DATA', Table({'IREFQ': unsigned, 'IMQ' : int, 'NOQ': int,
            'ITOQ CONC': Array(Sequence(int, float), as_list=tuple)
          }, numbering='IQ', row_condition=validate_row))
      ]
      super().__init__(name, members, has_hidden_members=True)

  def depends_on(self):
      return ['REFERENCE SYSTEM', 'MESH INFORMATION', 'TYPES']

  result_class = OccupationSection


section = OccupationSectionDefinition
