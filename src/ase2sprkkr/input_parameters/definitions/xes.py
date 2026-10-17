"""X-ray emission input, based on the common SPR-KKR 9.7 XRAYSPEC reader."""
import numpy as np
from ...common.grammar_types import Integer, Keyword, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ._file_manifests import XrayRateInputParametersDefinition
from ..input_parameters_definitions import InputValueDefinition as V
from ._sections import CONTROL, ENERGY, MODE, SITES, TASK, TAU


def _validate(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    if parameters.MODE.MODE() is not None:
        issues.append(DataValidityError("XES requires fully relativistic IREL=3."))
    shell = parameters.TASK.CL()
    if "SPDF".index(shell[1]) >= int(shell[0]):
        issues.append(DataValidityError("CL requires l < n."))
    if parameters.ENERGY.SPLITSS() or len(parameters.ENERGY.GRID()) != 1 or len(parameters.ENERGY.NE()) != 1:
        issues.append(DataValidityError("XES XRAYSPEC expects one TAU energy path."))
    for name in ("TAUCORE", "FRAMETET", "FRAMEPHI", "FRAMEGAM"):
        value = parameters.TASK[name]()
        if value is not None and not np.isfinite(value).all():
            issues.append(DataValidityError("XES widths and frame angles must be finite."))
    if parameters.TASK.TAUCORE() is not None and np.any(parameters.TASK.TAUCORE() < 0):
        issues.append(DataValidityError("TAUCORE must be nonnegative; zero delegates to the native width table."))
    if parameters.TASK.NONDIPE2() and parameters.TASK.NONDIPM1():
        issues.append(DataValidityWarning("NONDIPE2 takes precedence over NONDIPM1."))
    if parameters.TASK.ME() == "NAB" and parameters.TASK.DIPOL() == 0:
        issues.append(DataValidityWarning("ME=NAB forces DIPOL=1 in XRAYSPEC."))
    return issues or None


def input_parameters():
    energy = ENERGY(emax=(None, "XES forces the upper bound to EF", None), defaults={"GRID": 3, "NE": 300, "ImE": 0.01})
    energy.copy_member("EMAX", is_optional=True)
    energy.copy_member("NE", type=SetOf(Integer(min=2), min_length=1, max_length=2))
    task = TASK("XES").copy([
        V("IT", Integer(min=1), 1, is_always_added=False, info="Atom type, not site index."),
        V("CL", Keyword(*(f"{n}{s}" for n in range(1, 6) for s in "SPDF")), "2P", is_always_added=False),
        V("ME", Keyword("ADA", "NAB", "GRV"), "ADA", is_always_added=False),
        V("OUTPUT", Keyword("BARN", "MBARN", "SIGMA"), "BARN", is_always_added=False),
        V("POLOBS", Keyword("S", "C"), "S", is_always_added=False, info="Spherical (+,-,z) or Cartesian (x,y,z) polarization."),
        V("DIPOL", Integer(min=0, max=1), 0, is_always_added=False, info="0 includes non-dipole terms; NAB forces 1."),
        V("TAUCORE", SetOf(float, min_length=1, max_length=2), is_optional=True, info="Core widths in eV; last entry repeats. Read by XRAYSPEC, but no use in the RXES rate is established."),
        V("GFCONVRH", False), V("WRREGIRR", False),
        V("NONDIPE2", False), V("NONDIPM1", False),
        V("FRAMETET", float, is_optional=True, info="Frame Euler angle in degrees."),
        V("FRAMEPHI", float, is_optional=True, info="Frame Euler angle in degrees."),
        V("FRAMEGAM", float, is_optional=True, info="Frame Euler angle in degrees."),
    ])
    return XrayRateInputParametersDefinition("xes", [CONTROL("XES"), TAU, MODE, energy, task, SITES],
        executable="kkrgen", mpi=False, validators=_validate,
        info="X-ray emission through the common XRAYSPEC workflow.",
        description="Requires an occupied-state TAU directory with NE matching ENERGY.NE, or a workflow producing it. "
                    "Native XES forces EMAX=EF. RXES shares the calculation with RXAS; numerical emission "
                    "interpretation is unverified. MESPLIT is forced true and is not exposed as a tunable flag.")
