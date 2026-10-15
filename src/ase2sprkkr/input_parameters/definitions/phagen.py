"""PHAGEN input-parameter definition."""

from ._file_manifests import resolve_file
from .sections import CONTROL, TASK, TAU, ENERGY, SCF
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import TaskDefinition


class PhagenTaskDefinition(TaskDefinition):
    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        potfil = parameters.CONTROL.POTFIL()
        if potfil:
            files.extend(
                (
                    resolve_file(directory, f"{potfil}_ex_sprkkr"),
                    resolve_file(directory, f"{potfil}.only_coulomb"),
                )
            )
        return tuple(files)


def input_parameters_definition():
    schema = InputParametersDefinition(
        "phagen",
        [CONTROL("PHAGEN"), TASK("PHAGEN"), TAU, ENERGY(), SCF],
        info="PHAGEN",
        description="<TODO: change to something meaningfull>",
    )
    schema._task_definition = PhagenTaskDefinition(
        "PHAGEN",
        name="PHAGEN",
        executable="kkrscf",
        mpi=False,
    )
    return schema
