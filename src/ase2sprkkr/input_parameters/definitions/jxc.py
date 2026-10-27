"""JXC task definition."""

from ._file_manifests import resolve_file
from ._sections import CONTROL, TAU, ENERGY, SITES, STRCONST, MODE
from ..input_parameters_definitions import (
    InputParametersDefinition, InputSectionDefinition as Section,
    InputValueDefinition as V,
)
from ...common.generated_configuration_definitions import GeneratedValueDefinition
from ...common.grammar_types import flag, Integer, SetOf
from ...common.warnings import DataValidityError


class JxcInputParametersDefinition(InputParametersDefinition):
    task_selectors = ("JXC", "DMI")

    def task_definition_name(self, parameters):
        """The native flags identify DMI-only or JXC (also for combined runs)."""
        return _selected_task_name(parameters.TASK)

    def initialize_from_result(self, parameters, result):
        parameters.ENERGY.EMIN.set_default(result.last_iteration.energy.EMIN())

    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        dataset = parameters.CONTROL.DATASET()
        if parameters.TASK.JXC():
            if parameters.MODE.MODE() == "SP-SREL":
                files.extend((
                    resolve_file(directory, f"{dataset}_J_ij.dat"),
                    resolve_file(directory, f"{dataset}_Jij_*.agr"),
                ))
            else:
                files.extend((
                    resolve_file(directory, f"{dataset}_XCPLTEN_Jij.dat"),
                    resolve_file(directory, f"{dataset}_XCPLTEN_Jij_*.agr"),
                    resolve_file(directory, f"{dataset}_XCPLTEN_Dij.dat"),
                    resolve_file(directory, f"{dataset}_XCPLTEN_Dij_[xyz]_*.agr"),
                ))
        if parameters.TASK.DMI():
            files.append(resolve_file(directory, f"{dataset}_DMIVEC_Dij.dat"))
            files.append(resolve_file(directory, f"{dataset}_DMIVEC_Dij_[xyz]_*.agr"))
        return tuple(files)


def _requested_task_name(option):
    return option._get_root_container().requested_task_name


def _jxc_default(option):
    """Use the requested task for new inputs; infer native DMI-only on read."""
    requested = _requested_task_name(option)
    if requested is not None:
        return requested == "JXC"
    return not option._container.DMI()


def _dmi_default(option):
    return _requested_task_name(option) == "DMI"


def _selected_task_name(section):
    if section.JXC():
        return "JXC"
    if section.DMI():
        return "DMI"
    return "JXC"


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


def _validate_jxc(parameters, _values, why):
    if why == "set":
        return None

    issues = []
    jxc = parameters.TASK.JXC()
    dmi = parameters.TASK.DMI()
    mode = parameters.MODE.MODE()

    if not jxc and not dmi:
        issues.append(DataValidityError("At least one of TASK.JXC or TASK.DMI must be enabled."))
    if parameters.CONTROL.NONMAG():
        issues.append(DataValidityError("JXC/DMI requires a magnetic calculation."))
    if dmi and mode is not None:
        issues.append(DataValidityError("DMI requires magnetic FREL."))
    elif jxc and mode in ("NREL", "SREL"):
        issues.append(DataValidityError("JXC requires magnetic FREL or SP-SREL."))
    if parameters.ENERGY.GRID().tolist() != ["5"] or len(parameters.ENERGY.NE()) != 1 or parameters.ENERGY.SPLITSS():
        issues.append(DataValidityError("JXC/DMI requires one GRID=5 integration contour."))
    if parameters.TAU.BZINT() != "POINTS":
        issues.append(DataValidityError("The standard periodic JXC/DMI workflow requires BZINT=POINTS."))
    if parameters.TASK.SELECT() is not None and parameters.TASK.EXCLUDE() is not None:
        issues.append(DataValidityError("SELECT and EXCLUDE are mutually exclusive."))
    return issues or None


def input_parameters_definition():
    mmode = MODE["MODE"].copy(warning_condition=_jxc_mode_warning_condition)
    energy = ENERGY(defaults={"NE": [30], "ImE": 0.01})
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    tau = TAU.copy(
        add=[
            TAU["BZINT"].copy(default_value="POINTS"),
        ],
    )
    schema = JxcInputParametersDefinition(
        "jxc",
        [
            CONTROL("JXC"),
            tau,
            MODE.copy(add=[mmode]),
            STRCONST,
            energy,
            Section(
                "TASK",
                [
                    V("JXC", flag, default_value=_jxc_default, is_required=False,
                      info="Calculate exchange couplings and exchange tensor."),
                    V("DMI", flag, default_value=_dmi_default, is_required=False,
                      info="Calculate Dzyaloshinskii-Moriya interaction vectors."),
                    GeneratedValueDefinition("TASK", _selected_task_name),
                    V("HFF", False, info="Enable the hyperfine-field calculation."),
                    V("NSHLCLU", Integer(min=0), is_optional=True,
                      info="Shell count; overrides CLURAD."),
                    V("SELECT", SetOf(str, min_length=1, max_length=100),
                      is_optional=True, info="Included element symbols."),
                    V("EXCLUDE", SetOf(str, min_length=1, max_length=100),
                      is_optional=True, info="Excluded element symbols."),
                    V("QREPONLY", False, info="Treat only representative central sites."),
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
        info="JXC / DMI",
        validators=_validate_jxc,
        description="Periodic magnetic E-INTEGRAL JXC and DMI task with FREL "
                    "tensor couplings and DMIVECTOR, or SP-SREL scalar Jij.",
    )
    schema["CONTROL"]["NONMAG"].warning_condition = _jxc_nonmag_warning_condition
    return schema
