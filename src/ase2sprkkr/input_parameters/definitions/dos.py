"""DOS task definition."""

from .sections import CONTROL, TAU, ENERGY, TASK, SITES
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import TaskDefinition


def input_parameters_definition():
    return InputParametersDefinition(
        "dos",
        [
            CONTROL("DOS"),
            TAU,
            ENERGY(
                emax=(1.0, "value of the highest energy", None),
                defaults={"GRID": 3, "NE": 300, "ImE": 0.01},
            ),
            TASK("DOS"),
            SITES,
        ],
        info="DOS - The density of states computation",
    )


def task_definition():
    return TaskDefinition(
        input_parameters_definition,
        name="DOS",
        executable="kkrgen",
        mpi=True,
        output_file_templates=("{dataset}.dos",),
    )
