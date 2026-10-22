"""APS postprocessing of existing DOS (SPR-KKR 9.7 NRAPS)."""
import numpy as np
from ...common.grammar_types import Integer, Keyword
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters_definitions import InputValueDefinition as V
from ..task_definition import TaskDefinition
from ._file_manifests import NrSpectroscopyInputParametersDefinition
from .sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


def _validate_aps(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    shell = parameters.TASK.CL()
    if "SPDF".index(shell[1]) >= int(shell[0]):
        issues.append(DataValidityError("TASK.CL requires l < n."))
    if parameters.MODE.OP() in ("DMFT", "LDA+U"):
        issues.append(DataValidityError("APS does not support DMFT or LDA+U (GEN)."))
    if parameters.MODE.MODE() in (None, "SP-SREL"):
        issues.append(DataValidityWarning("APS forces IREL=min(1,IREL); use explicit SREL or NREL for clarity."))
    gamma = parameters.TASK.GAMMA()
    if gamma is not None:
        if not np.isfinite(gamma):
            issues.append(DataValidityError("TASK.GAMMA must be finite."))
        elif gamma <= 0:
            issues.append(DataValidityWarning("GAMMA<=0 disables DOS broadening."))
    if max(parameters.ENERGY.NE()) < 2:
        issues.append(DataValidityError("APS needs capacity for at least two DOS energies."))
    return issues or None


def input_parameters_definition():
    mode = MODE.copy()
    mode.copy_member("MODE", default_value="SREL")
    schema = NrSpectroscopyInputParametersDefinition(
        "aps",
        [CONTROL("DOS"), TAU, mode, STRCONST,
        ENERGY(emin=None, defaults={"GRID": 1, "NE": 300, "ImE": 0.0}),
        TASK("APS", [
            V("IT", Integer(min=1), 1, is_always_added=False),
            V("CL", Keyword(*(f"{n}{s}" for n in range(1, 6) for s in "SPDF")), "2P", is_always_added=False),
            V("CALCLS", False, info="Calculate (l,s)-resolved spectra."),
            V("OUTGNU", False, info="Write additional gnuplot data."),
            V("GAMMA", float, is_optional=True, info="Native parameter printed in eV; positive values enable DOS broadening."),
        ]), SITES],
        output_extension="aps",
        broadened_dos=True,
        validators=_validate_aps,
        info="Appearance-potential spectroscopy from existing DATASET.dos (ASA).",
        description="Requires old-format spin-polarized DOS on a uniform increasing mesh starting at EF. "
                    "ENERGY.NE supplies allocation capacity; actual energies come from DOS. "
                    "FULLPOT, DMFT and LDA+U are unsupported. Numerical output is unverified."
    )
    schema._task_definition = TaskDefinition("APS", executable="kkrgen", mpi=False, input_file_templates=("{dataset}_{adsi}.dos",))
    return schema
