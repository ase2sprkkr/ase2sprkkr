"""SCF input-parameter definition."""

from ._file_manifests import resolve_file
from .sections import CONTROL, TAU, ENERGY, SCF, SITES, STRCONST, CPA, MODE, TASK
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import TaskDefinition


class SCFInputParametersDefinition(InputParametersDefinition):
    """SCF schema with files depending on the concrete parameters."""

    def cache_files(self, parameters, directory):
        files = list(super().cache_files(parameters, directory))
        if parameters.SCF.FULLPOT():
            sfnfil = parameters.CONTROL.SFNFIL()
            if sfnfil:
                files.append(resolve_file(directory, sfnfil))
        return tuple(files)

    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        potfil = parameters.CONTROL.POTFIL()
        if potfil:
            files.append(resolve_file(directory, f"{potfil}_new"))
        return tuple(files)


def input_parameters_definition():
    schema = SCFInputParametersDefinition(
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
    )
    return schema
