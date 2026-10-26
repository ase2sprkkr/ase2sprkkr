"""Serial periodic ASA Landau susceptibility (SPR-KKR SUSC_LANDAU)."""
import numpy as np
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError
from ._sections import CONTROL, TAU, MODE, CPA, STRCONST, ENERGY, TASK, SITES
from ..input_parameters_definitions import InputParametersDefinition as InputParameters, InputValueDefinition as V, InputSectionDefinition as Section
from ..task_definition import TaskDefinition


def _validate_chilandau(p, _values, why):
    if why == "set":
        return None
    errors = []
    def error(message):
        errors.append(DataValidityError(message))
    if p.MODE.MODE() in ("NREL", "SREL", "SP-SREL"):
        error("This CHILANDAU definition covers FREL only.")
    if p.SCF.FULLPOT():
        error("Native CHILANDAU does not support FULLPOT.")
    if p.TAU.BZINT() != "POINTS":
        error("CHILANDAU requires BZINT=POINTS.")
    if p.ENERGY.GRID().tolist() != ["5"] or len(p.ENERGY.NE()) != 1 or p.ENERGY.SPLITSS():
        error("CHILANDAU requires one GRID=5 contour without SPLITSS; native code otherwise returns silently.")
    if not np.isfinite(p.TASK.CLURAD()) or p.TASK.CLURAD() < 0:
        error("CLURAD must be finite and nonnegative.")
    if p.TASK.SELECT() is not None and p.TASK.EXCLUDE() is not None:
        error("SELECT and EXCLUDE are mutually exclusive.")
    return errors or None


def input_parameters_definition():
    energy = ENERGY(defaults={"GRID": [5], "NE": [30], "ImE": 0.01})
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = InputParameters(
        "chilandau",
        [CONTROL("CHILANDAU"), TAU, MODE, CPA, STRCONST, energy,
            Section("SCF", [V("FULLPOT", False)], is_optional=True, is_expert=True),
            TASK("CHILANDAU").copy([
                V("CLURAD", 2.2, info="Neighbour radius in ALAT units; explicit convergence choice, not a native default."),
                V("NSHLCLU", Integer(min=0), is_optional=True, info="Shell count overrides CLURAD when supplied."),
                V("SELECT", SetOf(str, min_length=1, max_length=100), is_optional=True, info="Select neighbouring sites by element."),
                V("EXCLUDE", SetOf(str, min_length=1, max_length=100), is_optional=True, info="Exclude neighbouring sites by element."),
                V("QREPONLY", False, info="Representative central sites only; check multiplicities before enabling."),
                V("PRINT", Integer(min=0), is_optional=True, info="Native TASK print control."),
            ]), SITES],
        validators=_validate_chilandau,
        info="Landau orbital magnetic susceptibility",
        description="Serial periodic ASA FREL workflow; KMROT=0 required by native calculation. "
                    "Native MPI has allocation and reduction defects. Converge both energy mesh and neighbour radius. "
                    "Native output contains partial energy sums; use the last serial block."
    )
    schema._task_definition = TaskDefinition("CHILANDAU", experimental=True, executable="kkrchi", mpi=False)
    return schema
