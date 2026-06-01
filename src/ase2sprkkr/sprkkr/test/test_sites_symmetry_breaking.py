"""Tests for the symmetry-guard behaviour of Site proxies.

When a Site's SiteType is shared with other symmetry-equivalent sites,
any attempt to mutate a property through a proxy (or via a direct Site
setter) must raise a RuntimeError that directs the user to call
``site.break_symmetry()`` first.

After calling ``break_symmetry()`` the mutation must succeed and must
only affect the individual site — the other sites keep the original
shared SiteType.
"""

import pytest
import numpy as np

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

if True:
    from ..sprkkr_atoms import SPRKKRAtoms
    from ..radial_meshes import ExponentialMesh
    from ..radial import RadialPotential, RadialCharge
    from ..moments import Moments


def _two_shared_sites():
    """Return (atoms, s0, s1) where s0 and s1 share the same SiteType."""
    a = SPRKKRAtoms("NaCl")
    s0, s1 = a.sites[0], a.sites[1]
    s1.site_type = s0.site_type
    assert s0.site_type is s1.site_type
    return a, s0, s1


class TestSymmetryCheckOccupation(TestCase):

    def test_setter_raises_on_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        with pytest.raises(RuntimeError, match="break_symmetry"):
            s0.occupation = {"Fe": 1.0}

    def test_setter_works_after_break_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        s0.break_symmetry()
        s0.occupation = {"Fe": 1.0}
        self.assertEqual(s0.occupation.primary_symbol, "Fe")
        self.assertEqual(s1.occupation.primary_symbol, "Na")

    def test_setitem_raises_on_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        with pytest.raises(RuntimeError, match="break_symmetry"):
            s0.occupation["Na"] = 0.5

    def test_setitem_works_after_break_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        s0.break_symmetry()
        # add a second type so setitem can redistribute properly
        s0.occupation.add("Fe", 0.4)
        s0.occupation["Na"] = 0.7
        self.assertAlmostEqual(s0.occupation["Na"] + s0.occupation["Fe"], 1.0)
        self.assertAlmostEqual(s1.occupation["Na"], 1.0)

    def test_add_raises_on_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        with pytest.raises(RuntimeError, match="break_symmetry"):
            s0.occupation.add("Fe", 0.3)

    def test_add_works_after_break_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        s0.break_symmetry()
        s0.occupation.add("Fe", 0.3)
        self.assertEqual(len(s0.occupation), 2)
        self.assertEqual(len(s1.occupation), 1)

    def test_set_method_raises_on_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        with pytest.raises(RuntimeError, match="break_symmetry"):
            s0.occupation.set({"Fe": 0.6, "Ni": 0.4})

    def test_set_method_works_after_break_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        s0.break_symmetry()
        s0.occupation.set({"Fe": 0.6, "Ni": 0.4})
        self.assertEqual(s0.occupation.primary_symbol, "Fe")
        self.assertEqual(s1.occupation.primary_symbol, "Na")

    def test_delitem_raises_on_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        # set up two-element occupation on a solo site first, then re-share
        s0.break_symmetry()
        s0.occupation.set({"Na": 0.6, "Fe": 0.4})
        s1.site_type = s0.site_type
        with pytest.raises(RuntimeError, match="break_symmetry"):
            del s0.occupation["Fe"]

    def test_delitem_works_after_break_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        s0.break_symmetry()
        s0.occupation.set({"Na": 0.6, "Fe": 0.4})
        s1.site_type = s0.site_type
        s0.break_symmetry()
        del s0.occupation["Fe"]
        # After deletion, Vc fills the remainder so "Fe" is gone but Vc is present
        with self.assertRaises(KeyError):
            _ = s0.occupation["Fe"]
        # s1 keeps the original Na+Fe occupation
        self.assertEqual(len(s1.occupation), 2)

    def test_clean_raises_on_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        s0.break_symmetry()
        s0.occupation.set({"Na": 1.0, "Fe": 0.0})
        s1.site_type = s0.site_type
        with pytest.raises(RuntimeError, match="break_symmetry"):
            s0.occupation.clean()

    def test_clean_works_after_break_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        s0.break_symmetry()
        # Na:1.0, Fe:0.0 — sum is 1.0 so no Vc is added; clean() will drop Fe
        s0.occupation.set({"Na": 1.0, "Fe": 0.0})
        s1.site_type = s0.site_type
        s0.break_symmetry()
        s0.occupation.clean()
        self.assertEqual(len(s0.occupation), 1)
        self.assertEqual(len(s1.occupation), 2)

    def test_replace_type_raises_on_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        with pytest.raises(RuntimeError, match="break_symmetry"):
            s0.occupation.replace_type("Na", "Fe")

    def test_replace_type_works_after_break_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        s0.break_symmetry()
        s0.occupation.replace_type("Na", "Fe")
        self.assertEqual(s0.occupation.primary_symbol, "Fe")
        self.assertEqual(s1.occupation.primary_symbol, "Na")

    def test_reads_dont_raise(self):
        a, s0, s1 = _two_shared_sites()
        _ = s0.occupation.primary_symbol
        _ = s0.occupation["Na"]
        _ = len(s0.occupation)
        _ = list(s0.occupation)


class TestSymmetryCheckPotential(TestCase):

    def _pot(self, mesh, factor=1.0):
        # RadialPotential expects a 2D value array: row 0 = vt, row 1 = bt
        value = np.zeros((2, len(mesh.coors)))
        value[0] = mesh.coors * factor
        value[1] = mesh.coors * factor * 0.5
        return RadialPotential(value, mesh, 11)

    def test_setter_raises_on_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        m = ExponentialMesh(1.0, 0.005, None, None, 200, None)
        with pytest.raises(RuntimeError, match="break_symmetry"):
            s0.potential = self._pot(m)

    def test_setter_works_after_break_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        m = ExponentialMesh(1.0, 0.005, None, None, 200, None)
        p1 = self._pot(m, 1.0)
        p2 = self._pot(m, 2.0)
        s0.break_symmetry()
        s0.potential = p1
        s1.site_type = s0.site_type
        s0.break_symmetry()
        s0.potential = p2
        # potential(r) returns [vt(r), bt(r)]; check vt component at r=1.0
        self.assertAlmostEqual(s0.potential(1.0)[0], 2.0)
        self.assertAlmostEqual(s1.potential(1.0)[0], 1.0)

    def test_inplace_mutation_raises_on_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        m = ExponentialMesh(1.0, 0.005, None, None, 200, None)
        s0.break_symmetry()
        s0.potential = self._pot(m, 1.0)
        s1.site_type = s0.site_type
        original_bt = s0.potential.bt.copy()
        with pytest.raises(RuntimeError, match="break_symmetry"):
            s0.potential.bt = original_bt * 2

    def test_inplace_mutation_works_after_break_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        m = ExponentialMesh(1.0, 0.005, None, None, 200, None)
        s0.break_symmetry()
        s0.potential = self._pot(m, 1.0)
        s1.site_type = s0.site_type
        original_bt = s0.potential.bt.copy()
        s0.break_symmetry()
        s0.potential.bt = original_bt * 2
        np.testing.assert_allclose(s0.potential.bt, original_bt * 2)
        np.testing.assert_allclose(s1.potential.bt, original_bt)

    def test_reads_dont_raise(self):
        a, s0, s1 = _two_shared_sites()
        m = ExponentialMesh(1.0, 0.005, None, None, 200, None)
        s0.break_symmetry()
        s0.potential = self._pot(m)
        s1.site_type = s0.site_type
        _ = s0.potential
        _ = s0.potential.mesh


class TestSymmetryCheckMesh(TestCase):

    def test_setter_raises_on_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        new_mesh = ExponentialMesh(0.5, 0.005, None, None, 200, None)
        with pytest.raises(RuntimeError, match="break_symmetry"):
            s0.mesh = new_mesh

    def test_setter_works_after_break_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        original_mesh = s0.site_type.mesh
        new_mesh = ExponentialMesh(0.5, 0.005, None, None, 200, None)
        s0.break_symmetry()
        s0.mesh = new_mesh
        assert s0.mesh._get_target() is new_mesh
        assert s1.site_type.mesh is original_mesh

    def test_inplace_mutation_raises_on_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        with pytest.raises(RuntimeError, match="break_symmetry"):
            s0.mesh.r1 = s0.site_type.mesh.r1 * 2

    def test_inplace_mutation_works_after_break_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        original_r1 = s0.site_type.mesh.r1
        s0.break_symmetry()
        s0.mesh.r1 = original_r1 * 2
        self.assertAlmostEqual(s0.mesh.r1, original_r1 * 2)
        self.assertAlmostEqual(s1.site_type.mesh.r1, original_r1)

    def test_reads_dont_raise(self):
        a, s0, s1 = _two_shared_sites()
        _ = s0.mesh.r1
        _ = s0.mesh.dx


class TestSymmetryCheckMoments(TestCase):

    def _setup_with_moments(self):
        a, s0, s1 = _two_shared_sites()
        s0.break_symmetry()
        s0.moments = Moments(1.0, 2.0, 3.0, 4.0, 5.0)
        s1.site_type = s0.site_type
        return a, s0, s1

    def test_setter_raises_on_symmetry(self):
        a, s0, s1 = self._setup_with_moments()
        with pytest.raises(RuntimeError, match="break_symmetry"):
            s0.moments = Moments(9.0, 9.0, 9.0, 9.0, 9.0)

    def test_setter_works_after_break_symmetry(self):
        a, s0, s1 = self._setup_with_moments()
        s0.break_symmetry()
        s0.moments = Moments(9.0, 9.0, 9.0, 9.0, 9.0)
        self.assertAlmostEqual(s0.moments.smt, 9.0)
        self.assertAlmostEqual(s1.moments.smt, 3.0)

    def test_inplace_mutation_raises_on_symmetry(self):
        a, s0, s1 = self._setup_with_moments()
        with pytest.raises(RuntimeError, match="break_symmetry"):
            s0.moments.smt = 99.0

    def test_inplace_mutation_works_after_break_symmetry(self):
        a, s0, s1 = self._setup_with_moments()
        s0.break_symmetry()
        s0.moments.smt = 99.0
        self.assertAlmostEqual(s0.moments.smt, 99.0)
        self.assertAlmostEqual(s1.moments.smt, 3.0)

    def test_reads_dont_raise(self):
        a, s0, s1 = self._setup_with_moments()
        _ = s0.moments.smt


class TestSymmetryCheckReferenceSystem(TestCase):

    def test_setter_raises_on_symmetry(self):
        from ..reference_systems import ReferenceSystem
        a, s0, s1 = _two_shared_sites()
        original_vref = s0.site_type.reference_system.vref
        with pytest.raises(RuntimeError, match="break_symmetry"):
            s0.reference_system = ReferenceSystem(original_vref + 1.0, 0.0)

    def test_setter_works_after_break_symmetry(self):
        from ..reference_systems import ReferenceSystem
        a, s0, s1 = _two_shared_sites()
        original_vref = s0.site_type.reference_system.vref
        s0.break_symmetry()
        s0.reference_system = ReferenceSystem(original_vref + 1.0, 0.0)
        self.assertAlmostEqual(s0.reference_system.vref, original_vref + 1.0)
        self.assertAlmostEqual(s1.site_type.reference_system.vref, original_vref)

    def test_inplace_mutation_raises_on_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        original_vref = s0.site_type.reference_system.vref
        with pytest.raises(RuntimeError, match="break_symmetry"):
            s0.reference_system.vref = original_vref + 5.0

    def test_inplace_mutation_works_after_break_symmetry(self):
        a, s0, s1 = _two_shared_sites()
        original_vref = s0.site_type.reference_system.vref
        s0.break_symmetry()
        s0.reference_system.vref = original_vref + 5.0
        self.assertAlmostEqual(s0.reference_system.vref, original_vref + 5.0)
        self.assertAlmostEqual(s1.site_type.reference_system.vref, original_vref)

    def test_reads_dont_raise(self):
        a, s0, s1 = _two_shared_sites()
        _ = s0.reference_system.vref
