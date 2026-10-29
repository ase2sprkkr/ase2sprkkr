"""Stoner parameter input audited against SPR-KKR 9.7."""
import numpy as np
from ..task_definition import TaskDefinition
from ...common.warnings import DataValidityError
from ._file_manifests import TauCacheInputParametersDefinition
from ._sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


def _validate_stoner(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    energy = parameters.ENERGY
    if not np.array_equal(energy.GRID(), ["0"]) or not np.array_equal(energy.NE(), [1]):
        issues.append(DataValidityError("STONER requires one energy at EF: GRID=0, NE=1."))
    if energy.ImE() != 0 or energy.SPLITSS():
        issues.append(DataValidityError("STONER requires a real energy at EF and no SPLITSS."))
    for name in ("EMIN", "EMAX"):
        if energy[name]() is not None or energy[name + "EV"]() not in (None, 0):
            issues.append(DataValidityError("STONER energy bounds must delegate to EF or use a zero relative offset."))
    if parameters.MODE.MODE() is not None:
        issues.append(DataValidityError("STONER uses relativistic wave functions (WAVFUN_READ_REL); use FREL."))
    if parameters.SITES.NL() is not None and np.any(np.asarray(parameters.SITES.NL()) < 3):
        issues.append(DataValidityError("STONER accesses the d channel unconditionally; every NL must be at least 3."))
    return issues or None


def input_parameters_definition():
    energy = ENERGY(emin=(None, "Omitted bound delegates to EF", None),
                    emax=(None, "Omitted bound delegates to EF", None),
                    defaults={"GRID": 0, "NE": 1, "ImE": 0.0})
    for name in ("EMIN", "EMAX"):
        energy.copy_member(name, is_optional=True)
    schema = TauCacheInputParametersDefinition(
        "stoner",
        [CONTROL("STONER"), TAU, MODE, STRCONST, energy, TASK("STONER"), SITES],
        validators=_validate_stoner,
        info="Stoner parameter from the DOS and radial wave functions at EF.",
        description="No task-specific input options. Uses the VWN susceptibility kernel. "
                    "Native STONER enables reading a TAU directory and falls back to the normal KKRGEN calculation when no reusable TAU cache is available. "
                    "The final STONER branch is skipped by SPR-KKR in MPI mode. "
                    "Only the fully relativistic, single-energy workflow is exposed; "
                    "numerical results and FULLPOT have not been verified."
    )
    schema._task_definition = TaskDefinition("STONER", experimental=True, executable="kkrgen", mpi=False, output_files=(
            "{DATASET_ADSI}_*_g-spin.dat",
            "{DATASET_ADSI}_*_g-orb.dat",
            "{DATASET_ADSI}_*_VWN_A-xc.dat",
        ))
    return schema
