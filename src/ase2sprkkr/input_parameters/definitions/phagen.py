"""PHAGEN input-parameter definition."""

from ._file_manifests import resolve_file
from .sections import CONTROL, TASK, TAU, ENERGY, SCF
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import FileSpec, TaskDefinition


def _exchange_potential(parameters, directory):
    potfil = parameters.CONTROL.POTFIL()
    return resolve_file(directory, f"{potfil}_ex_sprkkr") if potfil else None


def _coulomb_potential(parameters, directory):
    potfil = parameters.CONTROL.POTFIL()
    return resolve_file(directory, f"{potfil}.only_coulomb") if potfil else None


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
        output_file_templates=(
            FileSpec(_exchange_potential),
            FileSpec(_coulomb_potential),
        ),
    )
    return schema
