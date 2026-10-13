from pathlib import Path

import pytest

from ase2sprkkr.potentials.potentials import Potential
from ase2sprkkr.sprkkr.calculator import SPRKKR


def _converged_potential():
    path = Path(__file__).parents[2] / "potentials" / "examples" / "FePt.new.pot"
    return Potential.from_file(path)


def test_save_input_restart_scf_preserves_canonical_potential(tmp_path):
    potential = _converged_potential()
    assert potential.SCF_INFO.SCFSTATUS() == "CONVERGED"

    _, _, potential_file, _ = SPRKKR().save_input(
        potential=potential,
        input_parameters="SCF",
        directory=tmp_path,
        restart_scf=True,
        empty_spheres=False,
    )

    assert potential.SCF_INFO.SCFSTATUS() == "CONVERGED"
    saved = Potential.from_file(potential_file)
    assert saved.SCF_INFO.SCFSTATUS() == "START"
    assert len(saved.POTENTIAL) == len(potential.POTENTIAL)
    assert len(saved.CHARGE) == len(potential.CHARGE)
    assert len(saved.MOMENTS) == len(potential.MOMENTS)


def test_restart_scf_is_rejected_for_non_scf_task(tmp_path):
    potential = _converged_potential()

    with pytest.raises(ValueError, match="only valid for an SCF"):
        SPRKKR().save_input(
            potential=potential,
            input_parameters="DOS",
            directory=tmp_path,
            restart_scf=True,
            empty_spheres=False,
        )
