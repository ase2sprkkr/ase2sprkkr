"""Experimental nuclear spin-spin coupling tensor, SPR-KKR 9.7."""
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError
from ..input_parameters_definitions import InputParametersDefinition as InputParameters, InputValueDefinition as V
from ..task_definition import TaskDefinition
from .sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


def _validate(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    if not parameters.CONTROL.NONMAG() or parameters.MODE.MODE() is not None:
        issues.append(DataValidityError("NUCSSCPL requires NONMAG and fully relativistic IREL=3."))
    if parameters.CONTROL.FINNUC():
        issues.append(DataValidityError("NUCSSCPL does not support FINNUC."))
    if parameters.ENERGY.GRID().tolist() != ["5"] or len(parameters.ENERGY.NE()) != 1 or parameters.ENERGY.SPLITSS():
        issues.append(DataValidityError("NUCSSCPL requires one GRID=5 energy path."))
    if parameters.TAU.BZINT() != "POINTS":
        issues.append(DataValidityError("NUCSSCPL off-diagonal TAU requires BZINT=POINTS."))
    if parameters.TASK.SELECT() is not None and parameters.TASK.EXCLUDE() is not None:
        issues.append(DataValidityError("SELECT and EXCLUDE are mutually exclusive."))
    return issues or None


def input_parameters_definition():
    control = CONTROL("NUCSSCPL").copy(add=[V("FINNUC", False, info="Finite nucleus; forbidden by this task.")])
    control.copy_member("NONMAG", default_value=True, is_always_added=True)
    tau = TAU.copy()
    tau.copy_member("BZINT", default_value="POINTS")
    energy = ENERGY(emax=(None, "GRID=5 uses the Fermi energy", None), defaults={"GRID": 5, "NE": 30, "ImE": 0.01})
    energy.copy_member("EMAX", is_optional=True)
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    task = TASK("NUCSSCPL").copy([
        V("CLURAD", 2.2, info="Cluster radius in ALAT units; explicit convenience default."),
        V("NSHLCLU", Integer(min=0), is_optional=True, info="Shell count; overrides CLURAD."),
        V("SELECT", SetOf(str, min_length=1, max_length=100), is_optional=True, info="Included element symbols."),
        V("EXCLUDE", SetOf(str, min_length=1, max_length=100), is_optional=True, info="Excluded element symbols."),
    ])
    schema = InputParameters(
        "nucsscpl",
        [control, tau, MODE, STRCONST, energy, task, SITES],
        validators=_validate,
        info="Experimental nuclear spin-spin coupling tensor (ASA).",
        description="Native 9.7 deallocates unallocated arrays on intermediate energies and clears "
                    "previous types' matrix elements. Fortran repair is needed before numerical use. "
                    "FULLPOT and KMROT!=0 are rejected by native code. Isotopes come from TABMUQIN; "
                    "the task does not read isotope or nuclear-moment input parameters."
    )
    schema._task_definition = TaskDefinition("NUCSSCPL", executable="kkrgen", mpi=False, output_file_templates=(
            "{dataset}_A_ij.dat",
            "RK.dat",
            "tab_aij.dat",
            "{dataset}_Aij_*.agr",
            "{dataset}_Rij_*.agr",
        ))
    return schema
