"""Experimental site-resolved back-scattering DOS, SPR-KKR 9.7."""
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError
from ..input_parameters_definitions import InputParametersDefinition as InputParameters, InputValueDefinition as V
from ..task_definition import TaskDefinition
from .sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


def _validate(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    if parameters.ENERGY.GRID().tolist() != ["3"] or len(parameters.ENERGY.NE()) != 1 or parameters.ENERGY.SPLITSS():
        issues.append(DataValidityError("BSDOSIJ requires a single GRID=3 energy path."))
    if parameters.MODE.MODE() is not None or parameters.CONTROL.NONMAG():
        issues.append(DataValidityError("BSDOSIJ requires magnetic FREL wave functions."))
    if parameters.TAU.BZINT() != "POINTS":
        issues.append(DataValidityError("BSDOSIJ requires BZINT=POINTS."))
    if parameters.TASK.SELECT() is not None and parameters.TASK.EXCLUDE() is not None:
        issues.append(DataValidityError("SELECT and EXCLUDE are mutually exclusive."))
    return issues or None


def input_parameters_definition():
    energy = ENERGY(emax=(None, "Omitted upper bound delegates to EF", None), defaults={"GRID": 3, "NE": 300, "ImE": 0.01})
    energy.copy_member("EMAX", is_optional=True)
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    tau = TAU.copy()
    tau.copy_member("BZINT", default_value="POINTS")
    task = TASK("BSDOSIJ").copy([
        V("CLURAD", 2.2, info="Cluster radius in ALAT units; explicit convenience default."),
        V("NSHLCLU", Integer(min=0), is_optional=True, info="Shell count; overrides CLURAD."),
        V("SELECT", SetOf(str, min_length=1, max_length=100), is_optional=True, info="Included element symbols."),
        V("EXCLUDE", SetOf(str, min_length=1, max_length=100), is_optional=True, info="Excluded element symbols."),
        V("QREPONLY", False, info="Treat only representative central sites."),
    ])
    schema = InputParameters(
        "bsdosij",
        [CONTROL("BSDOSIJ"), tau, MODE, STRCONST, energy, task, SITES],
        validators=_validate,
        info="Experimental site-resolved back-scattering DOS.",
        description="Native 9.7 component and output indexing appear faulty; numerical results are unverified. "
                    "Potential-dependent native checks also require KMROT=0, NKM=NKMMAX and representative type sites."
    )
    schema._task_definition = TaskDefinition("BSDOSIJ", executable="kkrgen", mpi=True, output_file_templates=(
            "{dataset}_{adsi}_XCPLTEN_Jij.dat",
            "{dataset}_{adsi}_XCPLTEN_Jij_*.agr",
            "{dataset}_{adsi}_XCPLTEN_Dij_[xyz]_*.agr",
        ))
    return schema
