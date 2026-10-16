"""SCF input-parameter definition."""

from .sections import CONTROL, TAU, ENERGY, SCF, SITES, STRCONST, CPA, MODE, TASK
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import FileSpec, TaskDefinition, TaskFiles


def _fullpot(parameters):
    return bool(parameters.SCF.FULLPOT())


def _sfn_file(parameters, _directory):
    return parameters.CONTROL.SFNFIL() or None


def _new_potential_file(parameters, _directory):
    potfil = parameters.CONTROL.POTFIL()
    return None if not potfil else f"{potfil}_new"


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
        files=TaskFiles(
            output=(FileSpec(_new_potential_file),),
            cache=(FileSpec(_sfn_file, condition=_fullpot),),
        ),
    )
    return schema
