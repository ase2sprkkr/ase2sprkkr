"""NRAES postprocessing of existing DOS (SPR-KKR 9.7 NRAES)."""
from ...common.grammar_types import Integer, Keyword
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters_definitions import InputParametersDefinition as InputParameters, InputValueDefinition as V
from ._sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


def _validate_nraes(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    shell = parameters.TASK.CL()
    if "SPDF".index(shell[1]) >= int(shell[0]):
        issues.append(DataValidityError("TASK.CL requires l < n."))
    if parameters.MODE.OP() in ("DMFT", "LDA+U"):
        issues.append(DataValidityError("NRAES does not support DMFT or LDA+U (GEN)."))
    if parameters.MODE.MODE() in (None, "SP-SREL"):
        issues.append(DataValidityWarning("NRAES forces IREL=min(1,IREL); use explicit SREL or NREL for clarity."))
    if max(parameters.ENERGY.NE()) < 2:
        issues.append(DataValidityError("NRAES needs capacity for at least two DOS energies."))
    return issues or None


def input_parameters():
    mode = MODE.copy()
    mode.copy_member("MODE", default_value="SREL")
    return InputParameters(
        "nraes", [CONTROL("NRAES"), TAU, mode, STRCONST,
        ENERGY(emin=None, defaults={"GRID": 1, "NE": 300, "ImE": 0.0}),
        TASK("NRAES", [
            V("IT", Integer(min=1), 1, is_always_added=False),
            V("CL", Keyword(*(f"{n}{s}" for n in range(1, 6) for s in "SPDF")), "2P", is_always_added=False),
            V("CALCLS", False, info="Calculate (l,s)-resolved spectra."),
            V("OUTGNU", False, info="Write additional gnuplot data."),
        ]), SITES], executable="kkrgen", mpi=False,
        input_files=("{dataset}.dos",), validators=_validate_nraes,
        info="Nonrelativistic Auger electron spectroscopy from existing DATASET.dos (ASA).",
        description="Requires old-format spin-polarized DOS on a uniform increasing mesh covering occupied states up to EF. "
                    "ENERGY.NE supplies allocation capacity; actual energies come from DOS. "
                    "FULLPOT, DMFT and LDA+U are unsupported. Numerical output is unverified.",
    )
