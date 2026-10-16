"""PHAGEN input-parameter definition."""

from .sections import CONTROL, TASK, TAU, ENERGY, SCF
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import FileSpec, TaskDefinition, TaskFiles


def _potential_output(parameters, suffix):
    potfil = parameters.CONTROL.POTFIL()
    return None if not potfil else f"{potfil}{suffix}"


def _extended_potential_file(parameters, _directory):
    return _potential_output(parameters, "_ex_sprkkr")


def _coulomb_potential_file(parameters, _directory):
    return _potential_output(parameters, ".only_coulomb")


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
                FileSpec(_extended_potential_file),
                FileSpec(_coulomb_potential_file),
            ),
        ),
    )
    return schema
