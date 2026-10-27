"""Tests for the symmetry-guard behaviour of Site proxies.

When a Site's SiteType is shared with other symmetry-equivalent sites,
any attempt to mutate a property through a proxy (or via a direct Site
setter) must raise a RuntimeError that directs the user to call
``site.break_symmetry()`` first.

After calling ``break_symmetry()`` the mutation must succeed and must
only affect the individual site — the other sites keep the original
shared SiteType.
"""

import numpy as np
import pytest

from ..moments import Moments
from ..radial import RadialPotential
from ..radial_meshes import ExponentialMesh
from ..sprkkr_atoms import SPRKKRAtoms


@pytest.fixture(scope="module")
def pristine_atoms():
    """Initialize the identical NaCl site state once per module."""
    atoms = SPRKKRAtoms("NaCl")
    _ = atoms.sites
    return atoms


@pytest.fixture
def shared_sites(pristine_atoms):
    """Copy the pristine state before each independent mutation scenario."""
    atoms = pristine_atoms.copy()
    site0, site1 = atoms.sites[0], atoms.sites[1]
    assert site0.site_type is not pristine_atoms.sites[0].site_type
    site1.site_type = site0.site_type
    assert site0.site_type is site1.site_type
    return atoms, site0, site1


def test_remesh_rejects_shared_site_then_changes_only_detached_site(shared_sites):
    _, site0, site1 = shared_sites
    original = site0.site_type.mesh
    mesh = ExponentialMesh(0.9, 0.005, None, None, 200, None)
    with pytest.raises(RuntimeError, match="break_symmetry"):
        site0.remesh(mesh)
    assert site0.site_type.mesh is site1.site_type.mesh is original
    site0.break_symmetry()
    site0.remesh(mesh)
    assert site0.site_type.mesh is mesh
    assert site1.site_type.mesh is original


class TestSymmetryCheckOccupation:
    def test_setter_rejects_shared_then_works_after_break_symmetry(self, shared_sites):
        _, site0, site1 = shared_sites
        assert site0.occupation.primary_symbol == "Na"
        assert site0.occupation["Na"] == pytest.approx(1.0)
        assert len(site0.occupation) == 1
        assert len(list(site0.occupation)) == 1
        with pytest.raises(RuntimeError, match="break_symmetry"):
            site0.occupation = {"Fe": 1.0}
        site0.break_symmetry()
        site0.occupation = {"Fe": 1.0}
        assert site0.occupation.primary_symbol == "Fe"
        assert site1.occupation.primary_symbol == "Na"

    def test_setitem_rejects_shared_then_works_after_break_symmetry(self, shared_sites):
        _, site0, site1 = shared_sites
        with pytest.raises(RuntimeError, match="break_symmetry"):
            site0.occupation["Na"] = 0.5
        site0.break_symmetry()
        site0.occupation.add("Fe", 0.4)
        site0.occupation["Na"] = 0.7
        assert site0.occupation["Na"] + site0.occupation["Fe"] == pytest.approx(1.0)
        assert site1.occupation["Na"] == pytest.approx(1.0)

    def test_add_rejects_shared_then_works_after_break_symmetry(self, shared_sites):
        _, site0, site1 = shared_sites
        with pytest.raises(RuntimeError, match="break_symmetry"):
            site0.occupation.add("Fe", 0.3)
        site0.break_symmetry()
        site0.occupation.add("Fe", 0.3)
        assert len(site0.occupation) == 2
        assert len(site1.occupation) == 1

    def test_set_method_rejects_shared_then_works_after_break_symmetry(self, shared_sites):
        _, site0, site1 = shared_sites
        with pytest.raises(RuntimeError, match="break_symmetry"):
            site0.occupation.set({"Fe": 0.6, "Ni": 0.4})
        site0.break_symmetry()
        site0.occupation.set({"Fe": 0.6, "Ni": 0.4})
        assert site0.occupation.primary_symbol == "Fe"
        assert site1.occupation.primary_symbol == "Na"

    def test_delitem_rejects_shared_then_works_after_break_symmetry(self, shared_sites):
        _, site0, site1 = shared_sites
        site0.break_symmetry()
        site0.occupation.set({"Na": 0.6, "Fe": 0.4})
        site1.site_type = site0.site_type
        with pytest.raises(RuntimeError, match="break_symmetry"):
            del site0.occupation["Fe"]
        site0.break_symmetry()
        del site0.occupation["Fe"]

        with pytest.raises(KeyError):
            _ = site0.occupation["Fe"]
        assert len(site1.occupation) == 2

    def test_clean_rejects_shared_then_works_after_break_symmetry(self, shared_sites):
        _, site0, site1 = shared_sites
        site0.break_symmetry()
        site0.occupation.set({"Na": 1.0, "Fe": 0.0})
        site1.site_type = site0.site_type
        with pytest.raises(RuntimeError, match="break_symmetry"):
            site0.occupation.clean()
        site0.break_symmetry()
        site0.occupation.clean()
        assert len(site0.occupation) == 1
        assert len(site1.occupation) == 2

    def test_replace_type_rejects_shared_then_works_after_break_symmetry(self, shared_sites):
        _, site0, site1 = shared_sites
        with pytest.raises(RuntimeError, match="break_symmetry"):
            site0.occupation.replace_type("Na", "Fe")
        site0.break_symmetry()
        site0.occupation.replace_type("Na", "Fe")
        assert site0.occupation.primary_symbol == "Fe"
        assert site1.occupation.primary_symbol == "Na"


class TestSymmetryCheckPotential:
    @staticmethod
    def _pot(mesh, factor=1.0):
        value = np.zeros((2, len(mesh.coors)))
        value[0] = mesh.coors * factor
        value[1] = mesh.coors * factor * 0.5
        return RadialPotential(value, mesh, 11)

    def test_setter_rejects_shared_then_works_after_break_symmetry(self, shared_sites):
        _, site0, site1 = shared_sites
        mesh = ExponentialMesh(1.0, 0.005, None, None, 200, None)
        potential1 = self._pot(mesh, 1.0)
        potential2 = self._pot(mesh, 2.0)
        with pytest.raises(RuntimeError, match="break_symmetry"):
            site0.potential = self._pot(mesh)
        site0.break_symmetry()
        site0.potential = potential1
        site1.site_type = site0.site_type
        site0.break_symmetry()
        site0.potential = potential2

        assert site0.potential(1.0)[0] == pytest.approx(2.0)
        assert site1.potential(1.0)[0] == pytest.approx(1.0)

    def test_inplace_mutation_rejects_shared_then_works_after_break_symmetry(self, shared_sites):
        _, site0, site1 = shared_sites
        mesh = ExponentialMesh(1.0, 0.005, None, None, 200, None)
        site0.break_symmetry()
        site0.potential = self._pot(mesh, 1.0)
        site1.site_type = site0.site_type
        assert site0.potential is not None
        assert site0.potential.mesh is not None
        original_bt = site0.potential.bt.copy()
        with pytest.raises(RuntimeError, match="break_symmetry"):
            site0.potential.bt = original_bt * 2
        site0.break_symmetry()
        site0.potential.bt = original_bt * 2
        np.testing.assert_allclose(site0.potential.bt, original_bt * 2)
        np.testing.assert_allclose(site1.potential.bt, original_bt)


class TestSymmetryCheckMesh:
    def test_setter_rejects_shared_then_works_after_break_symmetry(self, shared_sites):
        _, site0, site1 = shared_sites
        original_mesh = site0.site_type.mesh
        assert site0.mesh.r1 == original_mesh.r1
        assert site0.mesh.dx == original_mesh.dx
        new_mesh = ExponentialMesh(0.5, 0.005, None, None, 200, None)
        with pytest.raises(RuntimeError, match="break_symmetry"):
            site0.mesh = new_mesh
        site0.break_symmetry()
        site0.mesh = new_mesh
        assert site0.mesh._get_target() is new_mesh
        assert site1.site_type.mesh is original_mesh

    def test_inplace_mutation_rejects_shared_then_works_after_break_symmetry(self, shared_sites):
        _, site0, site1 = shared_sites
        original_r1 = site0.site_type.mesh.r1
        with pytest.raises(RuntimeError, match="break_symmetry"):
            site0.mesh.r1 = original_r1 * 2
        site0.break_symmetry()
        site0.mesh.r1 = original_r1 * 2
        assert site0.mesh.r1 == pytest.approx(original_r1 * 2)
        assert site1.site_type.mesh.r1 == pytest.approx(original_r1)


class TestSymmetryCheckMoments:
    @staticmethod
    def _setup_with_moments(shared_sites):
        atoms, site0, site1 = shared_sites
        site0.break_symmetry()
        site0.moments = Moments(1.0, 2.0, 3.0, 4.0, 5.0)
        site1.site_type = site0.site_type
        return atoms, site0, site1

    def test_setter_rejects_shared_then_works_after_break_symmetry(self, shared_sites):
        _, site0, site1 = self._setup_with_moments(shared_sites)
        assert site0.moments.smt == pytest.approx(3.0)
        with pytest.raises(RuntimeError, match="break_symmetry"):
            site0.moments = Moments(9.0, 9.0, 9.0, 9.0, 9.0)
        site0.break_symmetry()
        site0.moments = Moments(9.0, 9.0, 9.0, 9.0, 9.0)
        assert site0.moments.smt == pytest.approx(9.0)
        assert site1.moments.smt == pytest.approx(3.0)

    def test_inplace_mutation_rejects_shared_then_works_after_break_symmetry(self, shared_sites):
        _, site0, site1 = self._setup_with_moments(shared_sites)
        with pytest.raises(RuntimeError, match="break_symmetry"):
            site0.moments.smt = 99.0
        site0.break_symmetry()
        site0.moments.smt = 99.0
        assert site0.moments.smt == pytest.approx(99.0)
        assert site1.moments.smt == pytest.approx(3.0)


class TestSymmetryCheckReferenceSystem:
    def test_setter_rejects_shared_then_works_after_break_symmetry(self, shared_sites):
        from ..reference_systems import ReferenceSystem

        _, site0, site1 = shared_sites
        original_vref = site0.site_type.reference_system.vref
        assert site0.reference_system.vref == pytest.approx(original_vref)
        with pytest.raises(RuntimeError, match="break_symmetry"):
            site0.reference_system = ReferenceSystem(original_vref + 1.0, 0.0)
        site0.break_symmetry()
        site0.reference_system = ReferenceSystem(original_vref + 1.0, 0.0)
        assert site0.reference_system.vref == pytest.approx(original_vref + 1.0)
        assert site1.site_type.reference_system.vref == pytest.approx(original_vref)

    def test_inplace_mutation_rejects_shared_then_works_after_break_symmetry(self, shared_sites):
        _, site0, site1 = shared_sites
        original_vref = site0.site_type.reference_system.vref
        with pytest.raises(RuntimeError, match="break_symmetry"):
            site0.reference_system.vref = original_vref + 5.0
        site0.break_symmetry()
        site0.reference_system.vref = original_vref + 5.0
        assert site0.reference_system.vref == pytest.approx(original_vref + 5.0)
        assert site1.site_type.reference_system.vref == pytest.approx(original_vref)

