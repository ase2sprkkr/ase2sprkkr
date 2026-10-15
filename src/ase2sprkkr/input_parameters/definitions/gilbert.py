"""Gilbert input-parameter definition."""

from ...common.grammar_types import real, Keyword
from ._file_manifests import resolve_file
from .sections import CONTROL, TAU, ENERGY, SITES, TASK
from ..input_parameters_definitions import InputParametersDefinition, InputValueDefinition as V
from ..task_definition import TaskDefinition


class GilbertTaskDefinition(TaskDefinition):
    def input_files(self, parameters, directory):
        files = list(super().input_files(parameters, directory))
        if parameters.TASK.SETFLUCT() == "M_T" and parameters.TASK.FLUCTFIL():
            files.append(resolve_file(directory, parameters.TASK.FLUCTFIL()))
        return tuple(files)


def input_parameters_definition():
    schema = InputParametersDefinition(
        "gilbert",
        [
            CONTROL("GILBERT"),
            TAU,
            ENERGY(emin=None, emax=None, defaults={"GRID": 3, "NE": 1}),
            SITES.copy(defaults={"NL": 4}),
            TASK("Gilbert").copy(
                [
                    V("NTMP", int, is_required=False, info="Number of temperature points used for α(T)"),
                    V(
                        "SETFLUCT",
                        Keyword(
                            {
                                "MLIN": "Use linear temperature grid. Takes account only the electron scattering due to lattice vibrations.",
                                "M_T": "If NFTET = 1 and NFPHI = 1 but NVIBRA is bigger than 1, only lattice vibrations are taken into account. If NVIBRA = 1 but NFTET is bigger than 1 and NFPHI is bigger than 1, only spin fluctuations are taken into account.",
                            }
                        ),
                        default_value=None,
                        is_required=False,
                        info="Finite temperature calculation mode",
                    ),
                    V("NVIBRA", int, is_required=False, info="Number of directions for atomic displacements representing thermal lattice vibrations"),
                    V("TMPMIN", real, is_required=False, info="Lower limit of the temperature region (SETLFUNC=MLIN)."),
                    V("TMPMAX", real, is_required=False, info="Upper limit of the temperature region (SETLFUNC=MLIN)"),
                    V(
                        "FLUCTFIL",
                        str,
                        is_required=False,
                        info="Data file which contains temperature-dependent magnetization.",
                    ),
                    V("NFTET", int, is_required=False, info="Number of grid points specifying θ angle (SETLFUNC=M_T)"),
                    V("NFPHI", int, is_required=False, info="Number of grid points specifying φ angle (SETLFUNC=M_T)"),
                ]
            ),
        ],
        info="GILBERT - The Gilbert damping parameter calculation",
    )
    schema._task_definition = GilbertTaskDefinition(
        schema,
        name="GILBERT",
        executable="kkrchi",
    )
    return schema
