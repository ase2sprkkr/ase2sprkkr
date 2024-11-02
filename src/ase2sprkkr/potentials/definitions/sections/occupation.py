from ...potential_definitions import PotSectionDefinition, \
                                   PotValueDefinition
from ...potential_sections import PotentialSection

from ....common.grammar_types import Table, unsigned, Array, Sequence
import numpy as np
from ....sprkkr.sites import Site, Occupation
from ....sprkkr.spacegroup_info import SpacegroupInfo


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
      sg_info = SpacegroupInfo(atoms).recompute(init=True, atomic_kinds=indexes,
                                                update_info=False)
      unique = sg_info.dataset.equivalent_atoms if sg_info.dataset else \
               np.arange(len(indexes))

      tags = {}

      def site(i, d):
          ind = unique[i]
          if not ind in tags:
             occ = d['ITOQ CONC']
             mesh = read_io_data['meshes'][d['IMQ'] - 1]
             occ = Occupation({(read_io_data['types'][it - 1], oc) for it,oc in occ}, mesh=mesh)
             site = Site.create(atoms = atoms,
                         occupation = occ,
                         reference_system = read_io_data['reference_systems'][d['IREFQ'] - 1],
                         mesh = mesh)
             tags[ind] = site
             return site
          else:
             return Site(tags[ind].site_type)
      sites = np.array([ site(i,d) for i,d in enumerate(data()) ])
      atoms.set_sites(sites, sg_info)


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
