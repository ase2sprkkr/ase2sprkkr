import os

import ase.io
import numpy as np
import pytest

from ..calculator import SPRKKR
from ..radial import RadialCharge, RadialPotential
from ..radial_meshes import ExponentialMesh
from ..sites import SiteType
from ..sprkkr_atoms import SPRKKRAtoms


class TestSites:
    def test_adding(self):
        atoms = SPRKKRAtoms(symbols="NK", pbc=True)
        atoms.cell = [[1, 0, 0], [0, 1, 1], [0, 0, 1]]
        assert atoms.sites[0].site_type != atoms.sites[1].site_type
        out = atoms + atoms
        assert len(out) == 4
        assert id(out.sites[0]) != id(atoms.sites[0])
        assert id(out.sites[0].site_type) != id(atoms.sites[0].site_type)
        half = out[:2]
        assert id(out.sites[0]) != id(half.sites[0])
        assert id(out.sites[0].site_type) != id(half.sites[0].site_type)

    def test_potential(self):
        atoms = SPRKKRAtoms("NaCl")
        mesh1 = ExponentialMesh(1.0, 0.005, None, None, 200, None)
        mesh2 = ExponentialMesh(0.9, 0.005, None, None, 200, None)
        site0 = atoms.sites[0]
        site1 = atoms.sites[1]

        site0.potential = RadialPotential(mesh1.coors, mesh1, 2)
        site1.potential = site0.potential.copy()
        assert site0.potential(2.0) == 2.0
        assert site0.potential(1.0) == 1.0
        assert site0.potential(1.5) == 1.5
        assert site1.potential(2.0) == 2.0

        site1.potential = site0.potential.copy()
        assert site0.mesh != site1.mesh
        site1.remesh(site0.mesh)
        assert site0.mesh == site1.mesh

        site0.potential = potential2 = RadialPotential(mesh2.coors * 2, mesh2, 2)
        assert site0.potential(2.0) == 4.0
        assert site0.mesh == site1.mesh

        site0.mesh = mesh1
        assert site0.mesh == mesh1
        assert site0.potential.mesh == mesh1
        assert site0.potential.mesh == site0.mesh
        assert site0.potential(2.0) == 4.0
        assert site0.potential(1.0) == 2.0
        assert site0.potential(1.5) == 3.0

        assert site0.mesh != site1.mesh
        assert site0.potential.mesh != site1.mesh
        site1.remesh(site0.mesh)
        assert site0.mesh == site1.mesh
        assert site0.potential.mesh == site1.mesh
        site0.charge = mesh1.coors * 3
        assert site0.mesh == site1.mesh
        assert site0.charge(2.0) == 6.0
        assert site0.charge(1.0) == 3.0
        site0.potential = potential2
        assert site0.mesh != mesh2
        assert site0.potential.mesh != mesh2
        assert site0.charge.mesh != mesh2
        assert site0.charge(2.0) == 6.0
        assert site0.charge(1.0) == 3.0
        site0.charge = RadialCharge(site1.potential._value, mesh1)
        assert site0.mesh == mesh1
        assert site0.potential.mesh == mesh1
        assert site0.charge.mesh == mesh1

        site0.potential = [mesh1.coors, mesh1.coors * 2]
        np.testing.assert_allclose(site0.potential.interpolate(1.5), [1.5, 3.0])

        site0.potential = [1 / mesh1.coors, 1 / mesh1.coors * 2]
        np.testing.assert_allclose(site0.potential.interpolate(1.5), [1 / 1.5, 2 / 1.5])

    def test_vacuum(self):
        atoms = SPRKKRAtoms("NaCl")
        assert not SiteType(atoms, {"Na": 1.0}).is_vacuum()
        assert not SiteType(atoms, {"Na": 0.5, "Cl": 0.5}).is_vacuum()
        assert not SiteType(atoms, {"Na": 0.5, "Vc": 0.5}).is_vacuum()
        assert SiteType(atoms, {"Vc": 1.0}).is_vacuum()

        atoms = SPRKKRAtoms("NaX")
        assert not atoms.sites[0].is_vacuum()
        assert atoms.sites[1].is_vacuum()
        atoms.sites[1].atomic_types[0] = "Li"
        assert not atoms.sites[1].is_vacuum()

    def test_occupancy(self):
        atoms = SPRKKRAtoms("Na")
        site = atoms.sites[0]
        site.occupation = {"Na": 0.7}
        assert site.occupation[1] == pytest.approx(0.3)
        assert site.occupation.atomic_type(1).symbol == "Vc"
        site.occupation.atomic_type(1).symbol = "Fe"
        assert site.occupation.atomic_type(1).symbol == "Fe"
        assert site.occupation.atomic_type(1).n_electrons == 26

        atoms = SPRKKRAtoms("Na")
        site = atoms.sites[0]
        site.occupation = {"Na": 0.7, "Cl": None}
        assert site.occupation[1] == pytest.approx(0.3)
        assert site.occupation["Cl"] == pytest.approx(0.3)
        site.occupation["Cl"] = 0.6
        assert str(atoms.symbols) == "Cl"
        assert site.occupation[0] == pytest.approx(0.4)
        site.occupation["N"] = 0.5

        assert site.occupation[0] == pytest.approx(0.2)
        assert site.occupation[1] == pytest.approx(0.3)
        assert site.occupation[2] == pytest.approx(0.5)
        assert site.occupation["N"] == pytest.approx(0.5)
        assert str(atoms.symbols) == "N"
        site.occupation[0] = 0.0
        assert site.occupation[0] == 0.0
        assert site.occupation[1] == pytest.approx(0.3 * 5.0 / 4.0)
        assert site.occupation[2] == pytest.approx(0.5 * 5.0 / 4.0)
        assert len(site.occupation) == 3
        site.occupation.clean()
        assert len(site.occupation) == 2
        assert site.occupation[0] == pytest.approx(0.3 * 5.0 / 4.0)
        assert site.occupation[1] == pytest.approx(0.5 * 5.0 / 4.0)

    def test_occupancy_is_keyed_by_spacegroup_kind(self):
        atoms = SPRKKRAtoms("Na2")
        sites = atoms.sites
        sites[1].site_type = sites[0].site_type
        atoms.set_array("spacegroup_kinds", np.array([7, 7], dtype=int))
        atoms.info["occupancy"] = {}

        sites[0].site_type.occupation = {"Na": 0.75, "Cl": 0.25}

        assert set(atoms.info["occupancy"]) == {"7"}
        assert SiteType._ase_occupation(atoms, 0) == atoms.info["occupancy"]["7"]
        assert SiteType._ase_occupation(atoms, 1) == atoms.info["occupancy"]["7"]

    @pytest.mark.skip(reason="legacy slow integration test")
    def test_charge(self, calc_args):
        from ase.build import bulk

        atoms = bulk("Li")
        calculator = SPRKKR(atoms=atoms, **calc_args())
        calculator.input_parameters.SCF.MIX = 1.0

        for site in atoms.sites:
            site.mesh.jrws = 30
            site.mesh.r1 = 1e-6

        output = calculator.calculate(options={"NITER": 2}, atoms=atoms)
        atoms2 = output.atoms
        output = calculator.calculate(options={"NITER": 0}, atoms=atoms)
        atoms1 = output.atoms
        assert atoms2.sites[0].potential.vt != atoms1.sites[0].potential.vt
        output2 = calculator.calculate(atoms=atoms1, options={"NITER": 1})
        assert atoms2.sites[0].potential.vt == output2.atoms.sites[0].potential.vt
        assert atoms2.sites[0].potential.bt == output2.atoms.sites[0].potential.bt
        assert atoms2.sites[0].charge.raw_value == output2.atoms.sites[0].charge.raw_value
        atoms1.atoms.sites[0].potential.vt *= 2
        atoms1.atoms.sites[0].potential.bt *= 2
        output3 = calculator.calculate(atoms=atoms1, options={"NITER": 1})
        assert atoms2.sites[0].potential.vt != output3.atoms.sites[0].potential.vt

    def test_occupancy_and_spacegroup_kinds(self):
        atoms = ase.io.read(os.path.join(os.path.dirname(__file__), "test2.cif"))
        SPRKKRAtoms.promote_ase_atoms(atoms)
        assert str(atoms.symbols) == "X2I2FeX11"
        _ = atoms.sites
        assert str(atoms.symbols) == "X2I2FeX11"
