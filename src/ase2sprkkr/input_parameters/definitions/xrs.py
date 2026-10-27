"""X-ray resonant scattering input, based on the common SPR-KKR 9.7 XRAYSPEC reader."""
import numpy as np
from ...common.grammar_types import Integer, Keyword, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ._file_manifests import XrayRateInputParametersDefinition
from ..input_parameters_definitions import InputValueDefinition as V
from ..task_definition import TaskDefinition
from ._sections import CONTROL, ENERGY, MODE, SITES, TASK, TAU


def _validate(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    if parameters.MODE.MODE() is not None:
        issues.append(DataValidityError("XRS requires fully relativistic IREL=3."))
    shell = parameters.TASK.CL()
    if "SPDF".index(shell[1]) >= int(shell[0]):
        issues.append(DataValidityError("CL requires l < n."))
    if parameters.ENERGY.SPLITSS() or len(parameters.ENERGY.GRID()) != 1 or len(parameters.ENERGY.NE()) != 1:
        issues.append(DataValidityError("XRS XRAYSPEC expects one TAU energy path."))
    upper = parameters.ENERGY.EMAX()
    if not np.isfinite(parameters.TASK.EMAX3()) or (upper is not None and parameters.TASK.EMAX3() <= upper):
        issues.append(DataValidityError("EMAX3 must be finite and above ENERGY.EMAX."))
    for name in ("TAUCORE", "FRAMETET", "FRAMEPHI", "FRAMEGAM"):
        value = parameters.TASK[name]()
        if value is not None and not np.isfinite(value).all():
            issues.append(DataValidityError("XRS widths and frame angles must be finite."))
    if parameters.TASK.TAUCORE() is not None and np.any(parameters.TASK.TAUCORE() < 0):
        issues.append(DataValidityError("TAUCORE must be nonnegative; zero delegates to the native width table."))
    if parameters.TASK.NONDIPE2() and parameters.TASK.NONDIPM1():
        issues.append(DataValidityWarning("NONDIPE2 takes precedence over NONDIPM1."))
    if parameters.TASK.ME() == "NAB" and parameters.TASK.DIPOL() == 0:
        issues.append(DataValidityWarning("ME=NAB forces DIPOL=1 in XRAYSPEC."))
    return issues or None


def input_parameters_definition():
    energy = ENERGY(emin=None, emax=(4.0, "Explicit upper scan bound in Ry", None), defaults={"GRID": 6, "NE": 180, "ImE": 0.01})
    energy.copy_member("NE", type=SetOf(Integer(min=2), min_length=1, max_length=2))
    task = TASK("XRS", add=[
        V("IT", Integer(min=1), 1, is_always_added=False, info="Atom type, not site index."),
        V("CL", Keyword(*(f"{n}{s}" for n in range(1, 6) for s in "SPDF")), "2P", is_always_added=False),
        V("ME", Keyword("ADA", "NAB", "GRV"), "ADA", is_always_added=False),
        V("NE3", Integer(min=1), 40, is_always_added=False, info="Points on the third energy integration, GRID=7."),
        V("EMAX3", float, 8.0, is_always_added=False, info="Third integration upper bound in Ry."),
        V("TSELECT", Integer(min=0, max=2), 0, is_always_added=False, info="0 both terms; 1 resonant; 2 non-resonant."),
        V("POLOBS", Keyword("S", "C"), "S", is_always_added=False, info="Spherical (+,-,z) or Cartesian (x,y,z) polarization."),
        V("DIPOL", Integer(min=0, max=1), 0, is_always_added=False, info="0 includes non-dipole terms; NAB forces 1."),
        V("TAUCORE", SetOf(float, min_length=1, max_length=2), is_optional=True, info="Core widths in eV; last entry repeats. Read by XRAYSPEC, but no use in the RXRS rate is established."),
        V("GFCONVRH", False), V("WRREGIRR", False),
        V("NONDIPE2", False), V("NONDIPM1", False),
        V("FRAMETET", float, is_optional=True, info="Frame Euler angle in degrees."),
        V("FRAMEPHI", float, is_optional=True, info="Frame Euler angle in degrees."),
        V("FRAMEGAM", float, is_optional=True, info="Frame Euler angle in degrees."),
    ])
    schema = XrayRateInputParametersDefinition(
        "xrs",
        [CONTROL("XRS"), TAU, MODE, energy, task, SITES],
        validators=_validate,
        info="X-ray resonant scattering through the common XRAYSPEC workflow.",
        description="Three-integration resonant scattering tensor; requires NE matching the first TAU contour. "
                    "Native SIG is only partly initialized, so numerical results require Fortran repair. "
                    "OUTPUT is forced to SIGMA and MESPLIT to true. Dynamic workflow is not read for XRS."
    )
    schema._task_definition = TaskDefinition("XRS", experimental=True, executable="kkrgen", mpi=False)
    return schema
