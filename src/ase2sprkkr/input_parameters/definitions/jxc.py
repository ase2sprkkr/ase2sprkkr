"""JXC task definition."""

from ._file_manifests import resolve_file
from .sections import CONTROL, TAU, ENERGY, SITES, STRCONST, MODE
from ..input_parameters_definitions import (
    InputParametersDefinition, InputSectionDefinition as Section,
    InputValueDefinition as V,
)
from ...common.generated_configuration_definitions import GeneratedValueDefinition
from ...common.grammar_types import flag


class JxcInputParametersDefinition(InputParametersDefinition):
    task_selectors = ("JXC", "DMI")

    def task_definition_name(self, parameters):
        return _selected_task_name(parameters.TASK)

    def initialize_from_result(self, parameters, result):
        parameters.ENERGY.EMIN.set_default(result.last_iteration.energy.EMIN())

    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        dataset = parameters.CONTROL.DATASET()
        if parameters.TASK.JXC():
            if parameters.MODE.MODE() == "SP-SREL":
                files.append(resolve_file(directory, f"{dataset}*Jij*"))
            else:
                files.append(resolve_file(directory, f"{dataset}_XCPLTEN_Jij.dat"))
            files.append(resolve_file(directory, f"{dataset}*.agr"))
        if parameters.TASK.DMI():
            files.append(resolve_file(directory, f"{dataset}_DMIVEC_Dij.dat"))
            files.append(resolve_file(directory, f"{dataset}_DMIVEC_Dij_[xyz]_*.agr"))
        return tuple(files)


def _requested_task(option):
    return option._get_root_container().requested_task_name


def _jxc_default(option):
    requested = _requested_task(option)
    return requested == "JXC" or (
        requested is None and not option._container.DMI()
    )


def _dmi_default(option):
    return _requested_task(option) == "DMI"


def _selected_task_name(section):
    return "JXC" if section.JXC() else "DMI" if section.DMI() else "JXC"


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
            Section(
                "TASK",
                [
                    V("JXC", flag, default_value=_jxc_default, is_required=False,
                      info="Calculate exchange couplings and exchange tensor."),
                    V("DMI", flag, default_value=_dmi_default, is_required=False,
                      info="Calculate Dzyaloshinskii-Moriya interaction vectors."),
                    GeneratedValueDefinition("TASK", _selected_task_name),
                    V("HFF", False, info="Enable the hyperfine-field calculation."),
                    V(
                        "CLURAD",
                        2.2,
                        info="""The radius of a sphere restricting the cluster
        of atoms for which the exchange coupling
        constants will be calculated""",
                    ),
                ],
            ),
            SITES,
        ],
        info="JXC",
    )
    schema["CONTROL"]["NONMAG"].warning_condition = _jxc_nonmag_warning_condition
    return schema
