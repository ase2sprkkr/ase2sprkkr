"""JXC task definition."""

from ._file_manifests import resolve_file
from .sections import TASK, CONTROL, TAU, ENERGY, SITES, STRCONST, MODE
from ..input_parameters_definitions import InputParametersDefinition, InputValueDefinition as V
from ..task_definition import TaskDefinition
from ...common.grammar_types import flag


class JxcTaskDefinition(TaskDefinition):
    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        dataset = parameters.CONTROL.DATASET()
        if parameters.MODE.MODE() == "SP-SREL":
            files.append(resolve_file(directory, f"{dataset}*Jij*"))
        else:
            files.append(resolve_file(directory, f"{dataset}_XCPLTEN_Jij.dat"))
        files.append(resolve_file(directory, f"{dataset}*.agr"))
        return tuple(files)


def _jxc_mode_warning_condition(value):
    return (
        "JXC task does not support SREL (scalar relativity without spin) or "
        "NREL (no relativity at all) MODE. Please change SCF.MODE, or "
        "the computation will fail."
        if value in ("SREL", "NREL")
        else None
    )


def _jxc_nonmag_warning_condition(value):
    return (
        "JXC task does not support non-magnetic computation. Please disable "
        "CONTROL.NONMAG, or the computation will fail."
        if value
        else None
    )


def _prepare_from_result(parameters, result):
    parameters.set("EMIN", result.last_iteration.energy.EMIN())


def input_parameters_definition():
    mmode = MODE["MODE"].copy(warning_condition=_jxc_mode_warning_condition)
    schema = InputParametersDefinition(
        "jxc",
        [
            CONTROL("JXC"),
            TAU,
            MODE.copy(add=[mmode]),
            STRCONST,
            ENERGY(),
            TASK("JXC").copy(
                [
                    V(
                        "DMI",
                        flag,
                        default_value=False,
                        is_required=False,
                        info="Calculate Dzyaloshinskii-Moriya interaction vectors.",
                    ),
                    V(
                        "CLURAD",
                        2.2,
                        info="""The radius of a sphere restricting the cluster
        of atoms for which the exchange coupling
        constants will be calculated""",
                    ),
                ]
            ),
            SITES,
        ],
        info="JXC",
    )
    schema["CONTROL"]["NONMAG"].warning_condition = _jxc_nonmag_warning_condition
    schema._task_definition = JxcTaskDefinition(
        "JXC",
        name="JXC",
        executable="kkrgen",
        mpi=True,
        prepare_from_result=_prepare_from_result,
    )
    return schema
