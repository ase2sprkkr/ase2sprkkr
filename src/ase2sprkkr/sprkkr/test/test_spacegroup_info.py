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
