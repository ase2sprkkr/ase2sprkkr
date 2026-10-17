"""Shared JXC/DMI input-parameters definition."""

from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError
from ...common.generated_configuration_definitions import GeneratedValueDefinition

from ._file_manifests import resolve_file
from ._sections import TASK, CONTROL, TAU, ENERGY, SITES, STRCONST, MODE
from ..input_parameters_definitions import InputParametersDefinition, InputValueDefinition as V


class JxcInputParametersDefinition(InputParametersDefinition):
    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        dataset = parameters.CONTROL.DATASET()
        jxc = parameters.TASK.JXC()
        dmi = parameters.TASK.DMI()

        if jxc:
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
        if dmi:
            files.extend((
                resolve_file(directory, f"{dataset}_DMIVEC_Dij.dat"),
                resolve_file(directory, f"{dataset}_DMIVEC_Dij_[xyz]_*.agr"),
            ))
        return tuple(files)


def _task_name(_section):
    return "JXC"


def _requested_task_name(option):
    return getattr(option._container._get_root_container(), "_requested_task_name", "").upper()


def _jxc_default(option):
    return _requested_task_name(option) == "JXC"


def _dmi_default(option):
    return _requested_task_name(option) == "DMI"


def _mode_warning_condition(value):
    return (
        "JXC/DMI does not support SREL (scalar relativity without spin) or "
        "NREL (no relativity at all) MODE. Please change SCF.MODE, or "
        "the computation will fail."
        if value in ("SREL", "NREL")
        else None
    )


def _nonmag_warning_condition(value):
    return (
        "JXC/DMI does not support non-magnetic computation. Please disable "
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


def input_parameters():
    """JXC/DMI common E-INTEGRAL workflow."""
    mmode = MODE["MODE"].copy(warning_condition=_mode_warning_condition)
    energy = ENERGY(defaults={"NE": [30], "ImE": 0.01})
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    tau = TAU.copy()
    tau.copy_member("BZINT", default_value="POINTS")

    task = [
        V("JXC", bool, _jxc_default, info="Calculate exchange couplings/tensor."),
        V("DMI", bool, _dmi_default, info="Additionally or independently run the DMIVECTOR calculation."),
        GeneratedValueDefinition("TASK", _task_name),
        V("CLURAD", 2.2, info="Radius restricting the exchange-coupling cluster in ALAT units."),
        V("NSHLCLU", Integer(min=0), is_optional=True, info="Shell count; overrides CLURAD."),
        V("SELECT", SetOf(str, min_length=1, max_length=100), is_optional=True, info="Included element symbols."),
        V("EXCLUDE", SetOf(str, min_length=1, max_length=100), is_optional=True, info="Excluded element symbols."),
        V("QREPONLY", False, info="Treat only representative central sites."),
    ]
    task = TASK(task)

    input_parameters = JxcInputParametersDefinition(
        "jxc",
        [CONTROL("JXC"), tau, MODE.copy(add=[mmode]), STRCONST, energy, task, SITES],
        executable="kkrgen", mpi=True, info="JXC / DMI", validators=_validate_jxc,
        description="Standard periodic magnetic E-INTEGRAL workflow. JXC selects exchange couplings: "
                    "FREL tensor or SP-SREL scalar Jij. DMI selects DMIVECTOR and therefore requires FREL; "
                    "both flags may be enabled together. Native code requires KMROT=0; FREL also checks "
                    "NKM=NKMMAX and representative type sites. Embedded-cluster, occupation-weighted and "
                    "multisite higher-order workflows are not audited here.",
    )
    input_parameters["CONTROL"]["NONMAG"].warning_condition = _nonmag_warning_condition
    return input_parameters


aliases = ["dmi"]
