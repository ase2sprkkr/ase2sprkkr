import io
import os

from ...common.temporary import temporary_option_value
from ..potentials import Potential


def _converged_potential():
    path = os.path.join(os.path.dirname(__file__), "..", "examples", "FePt.new.pot")
    return Potential.from_file(path)


def test_restart_scf_is_only_a_serialization_override():
    potential = _converged_potential()
    assert potential.SCF_INFO.SCFSTATUS() == "CONVERGED"

    potential_count = len(potential.POTENTIAL)
    charge_count = len(potential.CHARGE)
    moments_count = len(potential.MOMENTS)

    output = io.StringIO()
    potential.save_to_file(output, restart_scf=True)

    assert potential.SCF_INFO.SCFSTATUS() == "CONVERGED"
    assert len(potential.POTENTIAL) == potential_count
    assert len(potential.CHARGE) == charge_count
    assert len(potential.MOMENTS) == moments_count

    output.seek(0)
    saved = Potential.from_file(output)
    assert saved.SCF_INFO.SCFSTATUS() == "START"
    assert len(saved.POTENTIAL) == potential_count
    assert len(saved.CHARGE) == charge_count
    assert len(saved.MOMENTS) == moments_count


def test_normal_save_preserves_start_status():
    potential = _converged_potential()
    potential.SCF_INFO.SCFSTATUS = "START"

    output = io.StringIO()
    potential.save_to_file(output)

    assert potential.SCF_INFO.SCFSTATUS() == "START"
    output.seek(0)
    assert Potential.from_file(output).SCF_INFO.SCFSTATUS() == "START"


def test_set_from_atoms_does_not_infer_scf_status():
    attached = _converged_potential()
    atoms = attached.atoms
    attached.SCF_INFO.SCFSTATUS = "ITR"

    standalone = Potential(atoms=atoms, set_from_atoms=False)
    assert standalone.SCF_INFO.SCFSTATUS() == "START"

    standalone.set_from_atoms(atoms)

    assert attached.SCF_INFO.SCFSTATUS() == "ITR"
    assert standalone.SCF_INFO.SCFSTATUS() == "START"
    assert atoms.potential is attached


def test_temporary_option_value_restores_existing_result():
    potential = _converged_potential()
    status = potential.SCF_INFO.SCFSTATUS
    status.result = "ITR"

    with temporary_option_value(status, "START"):
        assert status() == "START"
        assert status.result == "START"

    assert status() == "ITR"
    assert status.result == "ITR"
    status.clear_result()
    assert status() == "CONVERGED"


def test_restart_scf_restores_existing_result_override():
    potential = _converged_potential()
    status = potential.SCF_INFO.SCFSTATUS
    status.result = "ITR"

    output = io.StringIO()
    potential.save_to_file(output, restart_scf=True)

    assert status() == "ITR"
    assert status.result == "ITR"
    output.seek(0)
    assert Potential.from_file(output).SCF_INFO.SCFSTATUS() == "START"


def test_restart_override_does_not_change_intermediate_2d_status():
    potential = _converged_potential()
    potential.SCF_INFO.SCFSTATUS = "ITR-L-BULK"

    output = io.StringIO()
    potential.save_to_file(output, restart_scf=True)

    assert potential.SCF_INFO.SCFSTATUS() == "ITR-L-BULK"
    output.seek(0)
    assert Potential.from_file(output).SCF_INFO.SCFSTATUS() == "START"

    normal = io.StringIO()
    potential.save_to_file(normal)
    assert potential.SCF_INFO.SCFSTATUS() == "ITR-L-BULK"
    normal.seek(0)
    assert Potential.from_file(normal).SCF_INFO.SCFSTATUS() == "ITR-L-BULK"
