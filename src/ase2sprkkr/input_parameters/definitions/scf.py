"""SCF input-parameter definition."""

from .sections import CONTROL, TAU, ENERGY, SCF, SITES, STRCONST, CPA, MODE, TASK
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import TaskDefinition, TaskFile


def _shape_function_file(parameters):
    if not parameters.SCF.FULLPOT():
        return None
    return parameters.CONTROL.SFNFIL()


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
        output_files=(
            TaskFile(
                "{CONTROL.POTFIL}_new",
                info="Updated self-consistent potential",
            ),
        ),
        cache_files=(
            TaskFile(
                _shape_function_file,
                info="Full-potential shape-function cache file",
            ),
        ),
    )
    return schema
