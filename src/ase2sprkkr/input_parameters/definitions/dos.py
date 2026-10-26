"""DOS input-parameter definition."""

import numpy as np
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ._sections import CONTROL, TAU, ENERGY, TASK, SITES, MODE, CPA, STRCONST
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import TaskDefinition, TaskFile


def _validate_dos(p, _values, why):
    if why == "set":
        return None
    issues = []
    for grid, ne in zip(p.ENERGY.GRID(), p.ENERGY.NE()):
        if grid in ("2", "3") and ne < 2:
            issues.append(DataValidityError("DOS GRID=2/3 requires at least two energy points."))
    if not np.isfinite(p.ENERGY.ImE()) or p.ENERGY.ImE() < 0:
        issues.append(DataValidityError("DOS ImE must be finite and nonnegative."))
    if "3" in p.ENERGY.GRID() and p.ENERGY.ImE() == 0:
        issues.append(DataValidityWarning("DOS ImE=0 puts GRID=3 on the real axis; check poles and convergence."))
    for name in ("EMIN", "EMAX", "EMINEV", "EMAXEV"):
        value = p.ENERGY[name]()
        if value is not None and not np.isfinite(value):
            issues.append(DataValidityError(f"DOS {name} must be finite."))
    return issues or None


def input_parameters_definition():
    energy = ENERGY(
        emax=(1.0, "value of the highest energy", None),
        defaults={"GRID": 3, "NE": 300, "ImE": 0.01},
    )
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = InputParametersDefinition(
        "dos",
        [
            CONTROL("DOS"),
            TAU,
            energy,
            MODE,
            CPA,
            STRCONST,
            TASK("DOS"),
            SITES,
        ],
        info="DOS - The density of states computation",
        validators=_validate_dos,
    )
    schema._task_definition = TaskDefinition(
        "DOS",
        name="DOS",
        executable="kkrgen",
        mpi=True,
        output_reader="DOS",
        output_files=(TaskFile("{dataset}.dos", type="DOS"),),
    )
    return schema
