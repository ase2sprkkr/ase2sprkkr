from unittest.mock import patch

from ase.build import bulk
import numpy as np

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ...sprkkr.sprkkr_atoms import SPRKKRAtoms, SpacegroupInfo  # NOQA: E402
from ...sprkkr.spacegroup_info import SpacegroupInfoBlock  # NOQA: E402
from ...sprkkr.atoms_region import AtomsRegion  # NOQA: E402


class Test(TestCase):
    def test_spacegroup_info(self):

        def ssert(mapping, out):
            self.assertEqual(out, np.asarray(mapping))

        atoms = bulk("NaCl", "rocksalt", a=5.64)
        sgi = SpacegroupInfo(atoms)
        ssert([0, 1], sgi.equivalent_sites)
        self.assertEqual(225, sgi.spacegroup_number())
        self.assertEqual(225, sgi.dataset.number)

        atoms = SPRKKRAtoms("ION")
        atoms.set_positions(np.arange(9).reshape((3, 3)))
        sgi = SpacegroupInfo(atoms)
        assert sgi.dataset is False
        ssert([0, 1, 2], sgi.equivalent_sites)

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
            with self.assertRaises(RuntimeError):
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
            with self.assertRaises(NotImplementedError):
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
