"""Nuclear spin relaxation, with initialization safeguards for SPR-KKR 9.7."""
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters_definitions import InputParametersDefinition as InputParameters, InputValueDefinition as V
from ..task_definition import TaskDefinition
from .sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


def _validate_t1(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    mode = parameters.MODE.MODE()
    if mode in ("NREL", "SREL"):
        issues.append(DataValidityError("T1 NREL/SREL calls NONRELT1, which stops as not available."))
    elif mode == "SP-SREL":
        issues.append(DataValidityWarning("T1 SP-SREL enters RELT1; its relativistic channel assumptions are unverified."))
    if parameters.MODE.OP() in ("DMFT", "LDA+U"):
        issues.append(DataValidityError("T1 does not support DMFT or LDA+U (GEN)."))
    if parameters.CONTROL.PRINT() < 1:
        issues.append(DataValidityError("T1 requires PRINT>=1 to initialize NN/MM in SPR-KKR 9.7 RELT1."))
    if not parameters.CONTROL.RUNELOOP():
        issues.append(DataValidityError("T1 requires CONTROL.RUNELOOP to initialize TAUT before RELT1."))
    if parameters.ENERGY.NE().tolist() != [1] or parameters.ENERGY.GRID().tolist() != ["3"]:
        issues.append(DataValidityWarning("T1 forces one GRID=3 energy point at EF; GRID/NE selections are overwritten."))
    nl = parameters.SITES.NL()
    if nl is not None and max(nl) > 5:
        issues.append(DataValidityError("RELT1 has only 10 kappa labels; explicit SITES.NL must be <=5."))
    return issues or None


def input_parameters_definition():
    control = CONTROL("T1").copy([V("RUNELOOP", True, is_required=True,
                                    info="Override T1's native disabled energy loop to initialize TAUT.")])
    control.copy_member("PRINT", default_value=1,
                        info="T1 needs >=1 to initialize matrix-print dimensions in SPR-KKR 9.7; RELT1 later resets it to 0.")
    energy = ENERGY(emin=None, defaults={"GRID": 3, "NE": 1, "ImE": 0.0})
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = InputParameters(
        "t1",
        [control, TAU, MODE, STRCONST, energy, TASK("T1"), SITES],
        validators=_validate_t1,
        info="Nuclear spin relaxation at EF (ASA, serial; SPR-KKR numerical caveats).",
        description="No task-specific parameters are read. The energy is fixed at EF+i*ImE. "
                    "Nuclear data come from internal tables; no temperature or isotope option is read. "
                    "NREL/SREL, FULLPOT, DMFT and LDA+U are unavailable. PRINT=1 and RUNELOOP "
                    "are explicitly written to avoid native initialization defects. RELT1 additionally "
                    "uses signed kappa values as loop bounds, so numerical relaxation rates remain unverified."
    )
    schema._task_definition = TaskDefinition("T1", executable="kkrgen", mpi=False)
    return schema
