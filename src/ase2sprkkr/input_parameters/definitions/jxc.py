"""JXC task definition."""

from ._file_manifests import resolve_file
from .sections import TASK, CONTROL, TAU, ENERGY, SITES, STRCONST, MODE
from ..input_parameters_definitions import InputParametersDefinition, InputValueDefinition as V
from ..task_definition import FileSpec, TaskDefinition
from ...common.grammar_types import flag


def _jxc_spsrel_file(parameters, directory):
    dataset = parameters.CONTROL.DATASET()
    return resolve_file(directory, f"{dataset}*Jij*")


def _jxc_tensor_file(parameters, directory):
    dataset = parameters.CONTROL.DATASET()
    return resolve_file(directory, f"{dataset}_XCPLTEN_Jij.dat")


def _jxc_agr_file(parameters, directory):
    dataset = parameters.CONTROL.DATASET()
    return resolve_file(directory, f"{dataset}*.agr")


def _is_spsrel(parameters):
    return parameters.MODE.MODE() == "SP-SREL"


def _is_not_spsrel(parameters):
    return not _is_spsrel(parameters)


class JxcInputParametersDefinition(InputParametersDefinition):
    def update_from_result(self, parameters, result):
        parameters.set("EMIN", result.last_iteration.energy.EMIN())


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


def input_parameters_definition():
    mmode = MODE["MODE"].copy(warning_condition=_jxc_mode_warning_condition)
    schema = JxcInputParametersDefinition(
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
    schema._task_definition = TaskDefinition(
        "JXC",
        name="JXC",
        executable="kkrgen",
        mpi=True,
        output_reader="JXC",
        output_file_templates=(
            FileSpec(_jxc_spsrel_file, condition=_is_spsrel),
            FileSpec(_jxc_tensor_file, condition=_is_not_spsrel),
            FileSpec(_jxc_agr_file),
        ),
    )
    return schema
