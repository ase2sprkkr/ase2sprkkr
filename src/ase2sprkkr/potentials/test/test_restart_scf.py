import io
import os

from ...common.temporary import temporary_option_value
from ..potentials import Potential


def _converged_potential():
    path = os.path.join(os.path.dirname(__file__), "..", "examples", "FePt.new.pot")
    return Potential.from_file(path)


def test_restart_and_normal_save_preserve_status_and_electronic_data():
    potential = _converged_potential()
    assert potential.SCF_INFO.SCFSTATUS() == "CONVERGED"

    counts = {
        name: len(getattr(potential, name))
        for name in ("POTENTIAL", "CHARGE", "MOMENTS")
    }

    def serialize(*, restart_scf):
        output = io.StringIO()
        potential.save_to_file(output, restart_scf=restart_scf)
        output.seek(0)
        saved = Potential.from_file(output)
        for name, expected in counts.items():
            assert len(getattr(potential, name)) == expected
            assert len(getattr(saved, name)) == expected
        return saved.SCF_INFO.SCFSTATUS()

    # restart_scf affects only the written file, not the converged input.
    assert serialize(restart_scf=True) == "START"
    assert potential.SCF_INFO.SCFSTATUS() == "CONVERGED"

    # A normal write respects an explicit START state.
    potential.SCF_INFO.SCFSTATUS = "START"
    assert serialize(restart_scf=False) == "START"
    assert potential.SCF_INFO.SCFSTATUS() == "START"

    # The same contract holds for a non-converged 2D intermediate state.
    potential.SCF_INFO.SCFSTATUS = "ITR-L-BULK"
    assert serialize(restart_scf=True) == "START"
    assert potential.SCF_INFO.SCFSTATUS() == "ITR-L-BULK"
    assert serialize(restart_scf=False) == "ITR-L-BULK"
    assert potential.SCF_INFO.SCFSTATUS() == "ITR-L-BULK"


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


def test_temporary_override_and_restart_restore_existing_result():
    potential = _converged_potential()
    status = potential.SCF_INFO.SCFSTATUS
    status.result = "ITR"

    with temporary_option_value(status, "START"):
        assert status() == "START"
        assert status.result == "START"

    assert status() == "ITR"
    assert status.result == "ITR"

    output = io.StringIO()
    potential.save_to_file(output, restart_scf=True)

    assert status() == "ITR"
    assert status.result == "ITR"
    output.seek(0)
    assert Potential.from_file(output).SCF_INFO.SCFSTATUS() == "START"

    status.clear_result()
    assert status() == "CONVERGED"
