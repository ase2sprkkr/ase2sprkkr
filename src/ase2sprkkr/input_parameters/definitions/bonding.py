"""Experimental molecular bonding analysis (SPR-KKR 9.7)."""
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError
from ..input_parameters_definitions import InputParametersDefinition as InputParameters
from ..task_definition import TaskDefinition
from .sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


def _validate(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    if parameters.TAU.BZINT() != "CLUSTER" or not parameters.TAU.MOL():
        issues.append(DataValidityError("BONDING requires the molecular cluster workflow: BZINT=CLUSTER and MOL."))
    if parameters.MODE.MODE() is not None:
        issues.append(DataValidityError("BONDING hardcodes REL>CLM transformations; use FREL."))
    if parameters.ENERGY.SPLITSS() or len(parameters.ENERGY.GRID()) != 1 or len(parameters.ENERGY.NE()) != 1:
        issues.append(DataValidityError("BONDING returns without results for more than one energy path."))
    return issues or None


def input_parameters_definition():
    tau = TAU.copy()
    tau.copy_member("BZINT", default_value="CLUSTER")
    tau.copy_member("MOL", default_value=True, is_always_added=True)
    tau.copy_member("NSHLCLU", default_value=0, is_always_added=True)
    energy = ENERGY(emax=(1.0, "Explicit upper scan bound in Ry", None), defaults={"GRID": 3, "NE": 300, "ImE": 0.01})
    energy.copy_member("NE", type=SetOf(Integer(min=2), min_length=1, max_length=2))
    schema = InputParameters(
        "bonding",
        [CONTROL("BONDING"), tau, MODE, STRCONST, energy, TASK("BONDING"), SITES],
        validators=_validate,
        info="Experimental bonding analysis of molecular site 1.",
        description="Native 9.7 TAU_DRIVE supplies KTAUIJ=0, while CLUBONDING reads scratch units "
                    "that require KTAUIJ>=2. The native workflow needs repair before numerical use. "
                    "CLUBONDING also returns silently for KMROT!=0 and writes to unit 22."
    )
    schema._task_definition = TaskDefinition("BONDING", executable="kkrgen", mpi=False, output_file_templates=("fort.22",))
    return schema
