import io
import os
import pickle
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
from ase.spacegroup import crystal

from ...sprkkr.calculator import SPRKKR
from ..potential_sections import PotentialSection
from ..potentials import Potential


EXAMPLES = Path(__file__).parent.parent / "examples"


def _calculation_args(tmp_path, *, options, save_input=False):
    merged_options = {"NKTAB": 5, "NE": 20, **options}
    args = {
        "directory": tmp_path,
        "input_file": "output_test_calc.inp",
        "output_file": "output_test_calc.out",
        "potential_file": "output_test_calc.pot",
        "options": merged_options,
        "empty_spheres": False,
    }
    if not save_input:
        args.update(print_output=False, mpi="auto")
    return args


class TestPotential:
    def test_grammar(self):
        grammar = Potential.potential_definition.grammar()

        def check(expression, prefix=""):
            if "\n" in expression.whiteChars:
                raise AssertionError(
                    "The potential grammar has newline in whitespace, which is illegal:\n"
                    f" {expression}: {prefix}:{expression.whiteChars}"
                )
            if hasattr(expression, "expr"):
                check(expression.expr, prefix + ".expr")
            if hasattr(expression, "exprs"):
                for index, child in enumerate(expression.exprs):
                    check(child, prefix + f".exprs[{index}]")

        check(grammar)

    def test_potential_round_trip(self):
        a = 5.64
        atoms = crystal(
            ["Na", "Cl"],
            [(0, 0, 0), (0.5, 0.5, 0.5)],
            spacegroup=225,
            cellpar=[a, a, a, 90, 90, 90],
        )

        potential = Potential.from_atoms(atoms)
        potential.atoms.symbols[0] = "K"
        output = io.StringIO()
        potential.save_to_file(output)
        output.seek(0)
        restored = Potential.from_file(output)

        assert restored.atoms is not None
        np.testing.assert_allclose(restored.atoms.cell, potential.atoms.cell)
        np.testing.assert_allclose(restored.atoms.positions, potential.atoms.positions)
        assert str(restored.atoms.symbols) == str(potential.atoms.symbols)

    def test_scf_convergence_status(self):
        potential = Potential.from_file(EXAMPLES / "FePt.new.pot")
        atoms = potential.atoms

        assert atoms.sprkkr_is_scf_converged()

        potential.SCF_INFO.SCFSTATUS = "START"
        assert atoms.sprkkr_is_scf_converged()

        copied = atoms.copy()
        assert not copied.has_potential()
        assert copied.sprkkr_is_scf_converged()
        assert not copied.has_potential()

        for status in ("ITR", "ITR-BULK", "ITR-L-BULK", "ITR-R-BULK", "ITR-I-ZONE"):
            potential.SCF_INFO.SCFSTATUS = status
            assert not atoms.sprkkr_is_scf_converged()

            saved = io.StringIO()
            potential.save_to_file(saved)
            assert potential.SCF_INFO.SCFSTATUS() == status
            saved.seek(0)
            assert Potential.from_file(saved).SCF_INFO.SCFSTATUS() == status

        potential.SCF_INFO.SCFSTATUS = "START"
        potential.save_to_file(io.StringIO())
        assert potential.SCF_INFO.SCFSTATUS() == "CONVERGED"

    @pytest.mark.slow
    def test_2d(self, tmp_path):
        potential = Potential.from_file(EXAMPLES / "GeTe.pot")
        if os.environ.get("DO_NOT_RUN_SPRKKR", "") == "":
            SPRKKR().calculate(
                potential=potential,
                **_calculation_args(
                    tmp_path,
                    options={"NITER": 1, "NKTAB": 5, "NE": 5, "TOL": 1e-2},
                ),
            )
        else:
            SPRKKR().save_input(
                potential=potential,
                **_calculation_args(tmp_path, options={"NITER": 1}, save_input=True),
            )

    @pytest.mark.slow
    def test_potential_data(self, tmp_path):
        path = EXAMPLES / "FePt.new.pot"
        potential = Potential.from_file(path)
        assert potential.atoms.sites[0].potential is not None
        assert potential.atoms.sites[0].charge is not None
        assert potential.atoms.sites[-1].potential is not None
        assert potential.atoms.sites[-1].charge is not None

        potential.atoms.sites[1].potential.bt = -potential.atoms.sites[1].potential.bt
        saved = io.StringIO()
        potential.save_to_file(saved)
        saved.seek(0)
        restored = Potential.from_file(saved)
        original = Potential.from_file(path)

        np.testing.assert_allclose(
            restored.atoms.sites[1].potential.bt,
            -original.atoms.sites[1].potential.bt,
        )
        assert not np.allclose(
            restored.atoms.sites[1].potential.bt,
            original.atoms.sites[1].potential.bt,
        )

        if os.environ.get("DO_NOT_RUN_SPRKKR", "") == "":
            SPRKKR().calculate(
                potential=original,
                **_calculation_args(
                    tmp_path,
                    options={"NITER": 1, "NKTAB": 5, "NE": 5},
                ),
            )

    def test_examples(self):
        potential_files = list(EXAMPLES.glob("*.pot"))
        assert len(potential_files) >= 2

        for path in potential_files:
            potential = Potential.from_file(path)
            assert isinstance(potential, Potential)
            assert potential.atoms.positions.shape[0] > 0

    def test_sequence(self):
        potential = Potential()
        option = potential.HEADER.FORMAT
        option_class = type(option)

        assert potential.HEADER.FORMAT.VERSION == 9
        potential.HEADER.FORMAT.VERSION = 5
        assert potential.HEADER.FORMAT.VERSION == 5
        potential.HEADER.FORMAT[0] = 8
        assert potential.HEADER.FORMAT.VERSION == 8
        assert type(potential.HEADER.FORMAT.DATE) is datetime
        assert type(potential.HEADER.FORMAT[1]) is datetime

        restored = pickle.loads(pickle.dumps(option))
        assert type(restored) is option_class
        assert restored.VERSION == 8
        assert type(restored.DATE) is datetime

    def test_reset(self):
        potential = Potential.from_file(EXAMPLES / "fp_new.pot")
        potential.atoms.sites[0].mesh.r1 = 1.0
        potential.MESH_INFORMATION.DATA()[0][0] = 2.0
        assert potential.SCF_INFO.SCFSTATUS() != "START"

        # reset() clears non-mandatory potential properties but keeps the
        # mandatory structural information.
        potential.reset()
        assert potential.atoms.sites[0].mesh.r1 == 2.0
        assert potential.MESH_INFORMATION.DATA()[0][0] == 2.0
        assert len(potential.CHARGE) == 6
        assert len(potential.POTENTIAL) == 6
        assert isinstance(potential.LATTICE, PotentialSection)
        assert potential.SCF_INFO.SCFSTATUS() == "START"

        potential.atoms.sprkkr_reset_potential()
        assert potential.MESH_INFORMATION.DATA()[0][0] == potential.atoms.sites[0].mesh.r1
        assert potential.atoms.sites[0].mesh.r1 != 1.0
        assert potential.MESH_INFORMATION.DATA()[0][0] != 2.0
        assert potential.SCF_INFO.SCFSTATUS() == "START"
        assert len(potential.POTENTIAL) == 0
        assert len(potential.CHARGE) == 0
        assert len(potential.MOMENTS) == 0
        assert not potential.atoms.sprkkr_is_scf_converged()
        assert all(
            atomic_type.potential is None
            and atomic_type.charge is None
            and atomic_type.moments is None
            for site in potential.atoms.sites
            for atomic_type in site.occupation.atomic_types()
        )

    def test_single_type(self):
        potential = Potential.from_file(EXAMPLES / "Au_single_type.pot")
        assert potential.atoms is not None
