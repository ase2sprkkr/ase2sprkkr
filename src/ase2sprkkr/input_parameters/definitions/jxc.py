"""JXC task definition."""

from .sections import TASK, CONTROL, TAU, ENERGY, SITES, STRCONST, MODE
from ..input_parameters_definitions import InputParametersDefinition, InputValueDefinition as V
from ..task_definition import FileSpec, TaskDefinition, TaskFiles
from ...common.grammar_types import flag


def _jij_file(parameters, _directory):
    dataset = parameters.CONTROL.DATASET()
    if parameters.MODE.MODE() == "SP-SREL":
        return f"{dataset}*Jij*"
    return f"{dataset}_XCPLTEN_Jij.dat"


def _agr_file(parameters, _directory):
    return f"{parameters.CONTROL.DATASET()}*.agr"


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


class JxcInputParametersDefinition(InputParametersDefinition):
    def update_from_result(self, parameters, result):
        parameters.set("EMIN", result.last_iteration.energy.EMIN())


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
        files=TaskFiles(output=(FileSpec(_jij_file), FileSpec(_agr_file))),
    )
    return schema
