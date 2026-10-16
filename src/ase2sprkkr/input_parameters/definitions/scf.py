"""SCF input-parameter definition."""

from ._file_manifests import resolve_file
from .sections import CONTROL, TAU, ENERGY, SCF, SITES, STRCONST, CPA, MODE, TASK
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import FileSpec, TaskDefinition


def _sfn_file(parameters, directory):
    sfnfil = parameters.CONTROL.SFNFIL()
    return resolve_file(directory, sfnfil) if sfnfil else None


def _uses_fullpot(parameters):
    return parameters.SCF.FULLPOT()


def _new_potential_file(parameters, directory):
    potfil = parameters.CONTROL.POTFIL()
    return resolve_file(directory, f"{potfil}_new") if potfil else None


def input_parameters_definition():
    schema = InputParametersDefinition(
        "scf",
        [CONTROL("SCF"), TASK("SCF"), TAU, ENERGY(), SCF, SITES, STRCONST, CPA, MODE],
        info="SCF - Calculate a self-consistent potential",
        description="",
    )
    schema._task_definition = TaskDefinition(
        "SCF",
        name="SCF",
        executable="kkrscf",
        mpi=True,
        output_reader="SCF",
        output_file_templates=(FileSpec(_new_potential_file),),
        cache_file_templates=(FileSpec(_sfn_file, condition=_uses_fullpot),),
    )
    return schema
