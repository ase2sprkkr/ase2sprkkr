import numpy as np

from ..sites import Site
from ..sprkkr_atoms import SPRKKRAtoms


class TestSPRKKRAtomsCopyOwnership:
    def test_copy_sites_and_copy_preserve_sharing_without_reusing_source_types(self):
        atoms = SPRKKRAtoms("Fe2", positions=[[0, 0, 0], [1, 0, 0]], symmetry=False)
        atoms.sites[1].site_type = atoms.sites[0].site_type
        original_type = atoms.sites[0].site_type
        target = SPRKKRAtoms("Fe2", positions=atoms.positions, symmetry=False)

        copied_sites = Site.copy_sites(atoms.sites, target)
        assert copied_sites[0].site_type is copied_sites[1].site_type
        assert copied_sites[0].site_type.atoms is target
        assert copied_sites[0].site_type is not original_type

        copied = atoms.copy()
        assert copied.sites[0].site_type is copied.sites[1].site_type
        assert copied.sites[0].site_type is not original_type
        assert copied.sites[0].site_type.atoms is copied


    def test_copy_owns_independent_site_state(self):
        atoms = SPRKKRAtoms("Fe", symmetry=False)
        site = atoms.sites[0]
        site_type = site.site_type
        atomic_type = site_type.atomic_type(0)

        atomic_type.potential = np.ones((2, len(site_type.mesh)))
        atomic_type.charge = np.ones(len(site_type.mesh))
        atomic_type.moments = (1.0, 2.0, 3.0, 4.0, 5.0)
        atomic_type.magnetisation[:] = [1.0, 2.0, 3.0]
        site.magnetisation[:] = [4.0, 5.0]

        copied = atoms.copy()
        copied_site = copied.sites[0]
        copied_site_type = copied_site.site_type
        copied_atomic_type = copied_site_type.atomic_type(0)

        assert copied_site_type.atoms is copied
        assert copied_site_type is not site_type
        assert copied_site_type.mesh is not site_type.mesh
        assert copied_atomic_type is not atomic_type
        assert copied_atomic_type.mesh is copied_site_type.mesh
        assert copied_atomic_type.potential is not atomic_type.potential
        assert copied_atomic_type.charge is not atomic_type.charge
        assert copied_atomic_type.moments is not atomic_type.moments
        assert copied_atomic_type.potential.mesh is copied_site_type.mesh
        assert copied_atomic_type.charge.mesh is copied_site_type.mesh
        assert copied_atomic_type.magnetisation is not atomic_type.magnetisation
        assert copied_site.magnetisation is not site.magnetisation

        copied_atomic_type.potential.raw_value[0, 0] = 10.0
        copied_atomic_type.charge.raw_value[0] = 20.0
        copied_atomic_type.moments.smt = 30.0
        copied_atomic_type.magnetisation[0] = 40.0
        copied_site.magnetisation[0] = 50.0

        assert atomic_type.potential.raw_value[0, 0] == 1.0
        assert atomic_type.charge.raw_value[0] == 1.0
        assert atomic_type.moments.smt == 3.0
        assert atomic_type.magnetisation[0] == 1.0
        assert site.magnetisation[0] == 4.0
