"""X-ray absorption input, audited against SPR-KKR 9.7 XRAYSPEC."""
import numpy as np
from ...common.grammar_types import Integer, Keyword, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ._file_manifests import XrayRateInputParametersDefinition
from .sections import CONTROL, TASK, TAU, ENERGY, SITES, MODE
from ..input_parameters_definitions import InputValueDefinition as V
from ..task_definition import TaskDefinition


def _validate_xas(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    if parameters.MODE.MODE() is not None:
        issues.append(DataValidityError("XAS requires fully relativistic IREL=3."))
    shell = parameters.TASK.CL()
    if "SPDF".index(shell[1]) >= int(shell[0]):
        issues.append(DataValidityError("CL requires l < n."))
    if parameters.ENERGY.SPLITSS() or len(parameters.ENERGY.GRID()) != 1 or len(parameters.ENERGY.NE()) != 1:
        issues.append(DataValidityError("XAS XRAYSPEC expects one TAU energy path."))
    for name in ("FRAMETET", "FRAMEPHI", "FRAMEGAM"):
        value = parameters.TASK[name]()
        if value is not None and not np.isfinite(value):
            issues.append(DataValidityError("XAS frame angles must be finite."))
    if parameters.TASK.NONDIPE2() and parameters.TASK.NONDIPM1():
        issues.append(DataValidityWarning("NONDIPE2 takes precedence over NONDIPM1."))
    if parameters.TASK.ME() == "NAB" and parameters.TASK.DIPOL() == 0:
        issues.append(DataValidityWarning("ME=NAB forces DIPOL=1 in XRAYSPEC."))
    return issues or None


def input_parameters_definition():
    energy = ENERGY(emax=(4.0, "Upper bound in Ry; native XAS default", None),
                    defaults={"ImE": 0.01, "GRID": 6, "NE": 180})
    energy.copy_member("NE", type=SetOf(Integer(min=2), min_length=1, max_length=2))
    schema = XrayRateInputParametersDefinition(
        "xas",
        [CONTROL("XAS"), TAU, MODE, energy,
        TASK("XAS").copy([
            V("IT", Integer(min=1), 1, info="Atom type index, not site index."),
            V("CL", Keyword(*(f"{n}{s}" for n in range(1, 6) for s in "SPDF")), "2P", info="Initial core shell."),
            V("OUTPUT", Keyword("MBARN", "BARN", "SIGMA"), is_optional=True, info="Default cross section in Mb; SIGMA uses 10^15/s."),
            V("ME", Keyword("ADA", "NAB", "GRV"), "ADA", is_always_added=False),
            V("POLOBS", Keyword("S", "C"), "S", is_always_added=False),
            V("DIPOL", Integer(min=0, max=1), 0, is_always_added=False),
            V("GFCONVRH", False), V("WRREGIRR", False),
            V("NONDIPE2", False), V("NONDIPM1", False),
            V("FRAMETET", float, 0.0, is_expert=True, is_optional=True, info="Frame Euler polar angle in degrees."),
            V("FRAMEPHI", float, 0.0, is_expert=True, is_optional=True, info="Frame Euler azimuth in degrees."),
            V("FRAMEGAM", float, 0.0, is_expert=True, is_optional=True, info="Third frame Euler angle in degrees."),
        ]), SITES],
        validators=_validate_xas,
        info="Calculates X-ray absorption spectra",
        description="Requires NE to match the first TAU contour and a valid core shell in the potential. "
                    "Frame angles rotate the polarization frame. Native MPI exits before XRAYSPEC. "
                    "MECHECK is a separate task, not an XAS option; MESPLIT is forced true. "
                    "The DYNAMICAL alternative workflow is not covered by this definition."
    )
    schema._task_definition = TaskDefinition("XAS", executable="kkrgen", mpi=False)
    return schema
