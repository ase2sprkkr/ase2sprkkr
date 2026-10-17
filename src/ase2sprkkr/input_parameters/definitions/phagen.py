"""PHAGEN input-parameter definition."""

from .sections import CONTROL, TASK, TAU, ENERGY, SCF
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import TaskDefinition, TaskFile


def input_parameters_definition():
    schema = InputParametersDefinition(
        "phagen",
        [CONTROL("PHAGEN"), TASK("PHAGEN"), TAU, ENERGY(), SCF],
        info="PHAGEN",
        description="<TODO: change to something meaningfull>",
    )
    schema._task_definition = TaskDefinition(
        "PHAGEN",
        name="PHAGEN",
        executable="kkrscf",
        mpi=False,
        output_files=(
            TaskFile(
                "{CONTROL.POTFIL}_ex_sprkkr",
                info="Exchange part of the generated potential",
            ),
            TaskFile(
                "{CONTROL.POTFIL}.only_coulomb",
                info="Coulomb-only generated potential",
            ),
        ),
    )
    return schema
