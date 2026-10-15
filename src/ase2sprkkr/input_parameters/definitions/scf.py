"""SCF task definition."""

import os

from .sections import CONTROL, TAU, ENERGY, SCF, SITES, STRCONST, CPA, MODE, TASK
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import TaskDefinition


class SCFTaskDefinition(TaskDefinition):
    """SCF task with files depending on the concrete parameters."""

    @staticmethod
    def _path(directory, filename):
        return filename if os.path.isabs(filename) else os.path.join(directory, filename)

    def cache_files(self, parameters, directory):
        files = list(super().cache_files(parameters, directory))
        if parameters.SCF.FULLPOT():
            sfnfil = parameters.CONTROL.SFNFIL()
            if sfnfil:
                files.append(self._path(directory, sfnfil))
        return tuple(files)

    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        potfil = parameters.CONTROL.POTFIL()
        if potfil:
            files.append(self._path(directory, f"{potfil}_new"))
        return tuple(files)


def input_parameters_definition():
    return InputParametersDefinition(
        "scf",
        [CONTROL("SCF"), TASK("SCF"), TAU, ENERGY(), SCF, SITES, STRCONST, CPA, MODE],
        info="SCF - Calculate a self-consistent potential",
        description="",
    )


def task_definition():
    return SCFTaskDefinition(
        input_parameters_definition,
        name="SCF",
        executable="kkrscf",
        mpi=True,
    )
