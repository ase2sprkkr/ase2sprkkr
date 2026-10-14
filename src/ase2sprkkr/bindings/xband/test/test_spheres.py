import os

import ase.build
import pytest

from ase2sprkkr import SPRKKR


@pytest.mark.slow
@pytest.mark.skipif(
    bool(os.environ.get("DO_NOT_RUN_SPRKKR", "")),
    reason="SPR-KKR execution is disabled",
)
def test_empty_spheres_are_added_during_calculation(tmp_path):
    atoms = ase.build.bulk("Cu")

    SPRKKR().calculate(
        atoms,
        directory=str(tmp_path),
        input_file="output_test_calc.inp",
        output_file="output_test_calc.out",
        potential_file="output_test_calc.pot",
        mpi="auto",
        empty_spheres={"min_radius": 0.25},
        options={"NITER": 1, "NE": 20, "NKTAB": 5},
    )

    assert len(atoms) == 4
