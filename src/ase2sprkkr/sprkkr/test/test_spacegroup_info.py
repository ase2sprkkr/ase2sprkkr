from unittest.mock import patch

import numpy as np
import pytest
from ase.build import bulk

from ..atoms_region import AtomsRegion
from ..spacegroup_info import SpacegroupInfoBlock
from ..sprkkr_atoms import SPRKKRAtoms, SpacegroupInfo


class TestSpacegroupInfo:
    def test_spacegroup_info(self):
        atoms = bulk("NaCl", "rocksalt", a=5.64)
        sgi = SpacegroupInfo(atoms)
        np.testing.assert_array_equal(sgi.equivalent_sites, [0, 1])
        assert sgi.spacegroup_number() == 225
        assert sgi.dataset.number == 225

        atoms = SPRKKRAtoms("ION")
        atoms.set_positions(np.arange(9).reshape((3, 3)))
        sgi = SpacegroupInfo(atoms)
        assert sgi.dataset is False
        np.testing.assert_array_equal(sgi.equivalent_sites, [0, 1, 2])

    def test_nested_block_updating(self):
        sgi = SpacegroupInfo(SPRKKRAtoms("H"))

        with patch.object(SpacegroupInfoBlock, "recompute") as recompute:
            with sgi.block_updating():
                block = sgi._block
                assert block.counter == 1

                with sgi.block_updating():
                    assert sgi._block is block
                    assert block.counter == 2

                assert sgi._block is block
                assert block.counter == 1

            assert sgi._block is None
            recompute.assert_called_once_with(sgi)

    def test_block_updating_is_exception_safe(self):
        sgi = SpacegroupInfo(SPRKKRAtoms("H"))

        with patch.object(SpacegroupInfoBlock, "recompute") as recompute:
            with pytest.raises(RuntimeError):
                with sgi.block_updating():
                    raise RuntimeError("boom")

            assert sgi._block is None
            recompute.assert_called_once_with(sgi)

            with sgi.block_updating():
                pass
            assert sgi._block is None

    def test_block_updating_accumulates_recompute_options(self):
        sgi = SpacegroupInfo(SPRKKRAtoms("H"))

        with patch.object(SpacegroupInfoBlock, "recompute"):
            with sgi.block_updating():
                sgi.recompute(
                    init=True,
                    copy=True,
                    consider_old=False,
                    precision=1e-7,
                    angular_precision=2e-7,
                )
                block = sgi._block
                assert block.init
                assert block.copy
                assert block.update_info
                assert not block.consider_old
                assert block.precision == 1e-7
                assert block.angular_precision == 2e-7
                assert block.do

    def test_block_updating_rejects_atomic_kinds_array(self):
        sgi = SpacegroupInfo(SPRKKRAtoms("H"))

        with sgi.block_updating():
            with pytest.raises(NotImplementedError):
                sgi.recompute(atomic_kinds=np.array([0]))

    def test_update_spacegroup_kinds_if_required_repairs_invalid_state_only(self):
        atoms = SPRKKRAtoms("NaCl", symmetry=False)
        _ = atoms.sites
        sgi = atoms.spacegroup_info

        assert sgi.check_spacegroup_kinds()
        kinds = atoms.arrays["spacegroup_kinds"].copy()

        with patch.object(atoms, "set_array", wraps=atoms.set_array) as set_array:
            sgi.update_spacegroup_kinds(if_required=True)
            set_array.assert_not_called()

        kind = str(int(kinds[0]))
        del atoms.info["occupancy"][kind]
        assert not sgi.check_spacegroup_kinds()

        sgi.update_spacegroup_kinds(if_required=True)

        assert sgi.check_spacegroup_kinds()
        np.testing.assert_array_equal(atoms.arrays["spacegroup_kinds"], kinds)
        assert kind in atoms.info["occupancy"]

    def test_site_type_change_invalidates_dataset_and_rebuilds_kinds(self):
        atoms = SPRKKRAtoms("NaCl", symmetry=False)
        sites = atoms.sites
        atoms.spacegroup_info._dataset = object()

        sites[1].site_type = sites[0].site_type

        assert atoms.spacegroup_info._dataset is None
        assert atoms.spacegroup_info.check_spacegroup_kinds()
        np.testing.assert_array_equal(atoms.arrays["spacegroup_kinds"], [0, 0])
        assert list(atoms.info["occupancy"]) == ["0"]

    def test_region_info_is_cached_and_dataset_is_refreshed(self):
        atoms = SPRKKRAtoms(
            "HHe",
            positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.5]],
            cell=np.eye(3),
            pbc=False,
            symmetry=False,
        )
        region = AtomsRegion("left", slice(0, 1), np.eye(3), [False, False, False], atoms=atoms)
        sgi = atoms.spacegroup_info

        info = sgi.info_for_region(region)
        assert info is sgi.info_for_region(region)
        assert info.for_atoms is region
        assert info.dataset is False

        info._dataset = object()
        sgi.recompute()

        assert sgi.info_for_region(region) is info
        assert info.dataset is False


    def test_partial_sites_inherit_existing_symmetric_site_type(self):
        atoms = SPRKKRAtoms(
            bulk("Fe", "bcc", a=2.87, cubic=True), symmetry=True,
        )
        old = atoms.sites
        assert old[0].site_type is old[1].site_type
        original_type = old[1].site_type
        old[1].magnetisation[:] = [0.42, -0.19]
        partial = old.copy()
        partial[0] = 0  # first representative is unknown
        atoms.set_array("sprkkr_sites_data", partial)

        atoms.spacegroup_info.recompute(consider_old=True, copy=True)
        restored = atoms.sites
        assert restored[0].site_type is restored[1].site_type
        assert restored[0].site_type is not original_type
        assert restored[0].occupation.to_dict(True) == {"Fe": 1.0}
        np.testing.assert_allclose(restored[0].magnetisation, [0.42, -0.19])
        np.testing.assert_allclose(restored[1].magnetisation, [0.42, -0.19])
        assert restored[1] is not old[1]

    def test_periodic_regions_stay_distinct_in_global_search(self):
        atoms = SPRKKRAtoms(
            bulk("Fe", "bcc", a=2.87, cubic=True), symmetry=True,
        )
        # The regions are processed first, then the whole structure. The
        # latter must not merge their separately assigned Sites.
        for i in range(2):
            region = AtomsRegion(
                f"region_{i}", slice(i, i + 1),
                atoms.cell.array, [True, True, True],
                atoms=atoms, add=False,
            )
            atoms.add_region(region)
        assert atoms.sites[0].site_type is not atoms.sites[1].site_type

    def test_explicit_kinds_can_merge_distinct_compositions(self):
        atoms = SPRKKRAtoms(
            bulk("NaCl", "rocksalt", a=5.64), symmetry=True,
        )
        atoms.spacegroup_info.recompute(
            atomic_kinds=np.array([2, 2]), symmetry=True,
        )
        assert atoms.sites[0].site_type is atoms.sites[1].site_type


    def test_initial_magmoms_keep_two_bcc_sites_inequivalent(self):
        atoms = SPRKKRAtoms(bulk("Fe", "bcc", a=2.87, cubic=True))
        atoms.set_initial_magnetic_moments([2.0, -2.0])
        assert atoms.sites[0].site_type is not atoms.sites[1].site_type
        assert len(np.unique(atoms.spacegroup_info.dataset.equivalent_atoms)) == 2

    def test_existing_site_magnetisation_survives_recompute(self):
        atoms = SPRKKRAtoms(
            bulk("Fe", "bcc", a=2.87, cubic=True), symmetry=True,
        )
        old_sites = atoms.sites
        assert old_sites[0].site_type is old_sites[1].site_type
        old_sites[0].magnetisation[:] = [0.6, 0.1]
        old_sites[1].magnetisation[:] = [-0.2, 0.8]

        atoms.spacegroup_info.recompute(consider_old=True, copy=True)
        new_sites = atoms.sites
        assert new_sites[0].site_type is not new_sites[1].site_type
        np.testing.assert_allclose(new_sites[0].magnetisation, [0.6, 0.1])
        np.testing.assert_allclose(new_sites[1].magnetisation, [-0.2, 0.8])
        assert new_sites[0] is not old_sites[0]
        assert new_sites[1] is not old_sites[1]

    def test_regions_separate_memberships_but_preserve_internal_equivalence(self):
        atoms = SPRKKRAtoms(
            bulk("Fe", "bcc", a=2.87, cubic=True).repeat((2, 1, 1)),
            symmetry=True,
        )
        for name, selector in (
            ("left", slice(0, 2)), ("right", slice(2, 4)),
        ):
            atoms.add_region(
                AtomsRegion(
                    name, selector, atoms.cell.array, [True, True, True],
                    atoms=atoms, add=False,
                )
            )
        sites = atoms.sites
        # Symmetry within a region must remain, between regions must not.
        assert sites[0].site_type is sites[1].site_type
        assert sites[2].site_type is sites[3].site_type
        assert sites[0].site_type is not sites[2].site_type
        assert atoms.spacegroup_info.check_spacegroup_kinds()

    def test_explicit_atomic_kinds_with_initial_magmoms_still_override(self):
        atoms = SPRKKRAtoms(bulk("Fe", "bcc", a=2.87, cubic=True))
        atoms.set_initial_magnetic_moments([1.0, -1.0])
        atoms.spacegroup_info.recompute(atomic_kinds=np.array([8, 8]))
        assert atoms.sites[0].site_type is atoms.sites[1].site_type


    def test_recompute_unregisters_replaced_sites(self):
        atoms = SPRKKRAtoms(
            bulk("Fe", "bcc", a=2.87, cubic=True), symmetry=False,
        )
        original = atoms.sites.copy()
        original_types = [site.site_type for site in original]
        atoms.spacegroup_info.recompute(consider_old=True, copy=False)
        assert all(new is not old for new, old in zip(atoms.sites, original))
        for site, site_type in zip(original, original_types):
            assert site not in site_type.sites
        for site in atoms.sites:
            assert site in site.site_type.sites
            assert len(site.site_type.sites) == 1
