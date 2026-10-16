"""PHAGEN input-parameter definition."""

from .sections import CONTROL, TASK, TAU, ENERGY, SCF
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import FileSpec, TaskDefinition, TaskFiles


def _potential_output(suffix):
    def resolve(parameters, _directory):
        potfil = parameters.CONTROL.POTFIL()
        return None if not potfil else f"{potfil}{suffix}"

    return resolve


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
        files=TaskFiles(
            output=(
                FileSpec(_potential_output("_ex_sprkkr")),
                FileSpec(_potential_output(".only_coulomb")),
            ),
        ),
    )
    return schema
