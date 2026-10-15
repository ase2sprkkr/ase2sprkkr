import os

import numpy as np
import pytest
from ase.atoms import Atoms

from ..calculator import SPRKKR


@pytest.mark.slow
def test_arpes_run(calc_args):
    if os.environ.get("DO_NOT_RUN_SPRKKR", ""):
        pytest.skip("SPR-KKR execution disabled by DO_NOT_RUN_SPRKKR")

    atoms = Atoms(
        symbols="Cu",
        cell=np.array([(1.0, 0, 0), (0, 1, 0), (0, 0, 1)]),
        pbc=True,
    )
    calculator = SPRKKR(atoms=atoms, **calc_args())
    calculator.input_parameters.SCF.NITER = 1
    calculator.set("NE", 10)
    calculator.set("NKTAB", 10)

    scf_output = calculator.calculate(**calc_args())
    arpes_calculator = scf_output.calculator
    arpes_calculator.input_parameters = "arpes"
    arpes_calculator.set("NE", 10)
    arpes_calculator.set("NKTAB", 10)

    arpes_output = arpes_calculator.calculate(**calc_args())
    assert isinstance(arpes_output.spc.ENERGY(), np.ndarray)
