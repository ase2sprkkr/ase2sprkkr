import pytest
from ase.build import bulk

from ...sprkkr.calculator import SPRKKR
from ...sprkkr.sprkkr_atoms import SPRKKRAtoms
from ..potentials import Potential


def test_occupancy_is_written_to_potential(tmp_path):
    atoms = SPRKKRAtoms.promote_ase_atoms(bulk("Cu", a=3.6))
    atoms.sites[0].occupation.set({"Cu": 0.75, "Ni": 0.25})

    options = {
        "SITES.NL": 3,
        "TAU.BZINT": "POINTS",
        "TAU.NKTAB3D": 500,
        "SCF.VXC": "VWN",
        "SCF.NITER": 200,
        "SCF.MIX": 0.1,
        "SCF.MIXOP": 0.1,
        "SCF.TOL": 1e-5,
    }

    calculator = SPRKKR(atoms=atoms)
    calculator.save_input(
        input_file="Cu.inp",
        potential_file="Cu.pot",
        options=options,
        directory=tmp_path,
    )

    assert (tmp_path / "Cu.inp").is_file()
    potential_file = tmp_path / "Cu.pot"
    assert potential_file.is_file()

    restored = Potential.from_file(potential_file).atoms.sites[0].occupation
    concentrations = {str(symbol): concentration for symbol, concentration in restored.items()}
    assert concentrations == pytest.approx({"Cu": 0.75, "Ni": 0.25})
