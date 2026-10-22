"""FSCAT input for the fixed backscattering geometries in SPR-KKR 9.7."""
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters_definitions import InputParametersDefinition as InputParameters
from ..task_definition import TaskDefinition
from .sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


def _validate_fscat(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    if parameters.MODE.OP() in ("DMFT", "LDA+U"):
        issues.append(DataValidityError("FSCAT does not support DMFT or LDA+U (GEN)."))
    if parameters.MODE.MODE() is not None:
        issues.append(DataValidityWarning(
            "FSCAT constructs relativistic kappa/spin channels; other MODE settings are unverified."))
    energy = parameters.ENERGY
    if energy.ImE() != 0:
        issues.append(DataValidityWarning("FSCAT resets ENERGY.ImE to zero."))
    if energy.SPLITSS():
        issues.append(DataValidityWarning("FSCAT disables ENERGY.SPLITSS."))
    if len(energy.NE()) > 1 or len(energy.GRID()) > 1:
        issues.append(DataValidityWarning("FSCAT uses only the first ENERGY path."))
    nl = parameters.SITES.NL()
    if nl is not None:
        if max(nl) < 2:
            issues.append(DataValidityError(
                "FSCAT needs NL>=2: CALC_CHPLM writes l=1 entries even for an l_max=0 allocation."))
        if len(set(nl)) > 1:
            issues.append(DataValidityWarning("FSCAT overwrites local cutoffs with the global maximum NL."))
    return issues or None


def input_parameters_definition():
    energy = ENERGY(emax=(None, "Upper energy in Ry; omitted means the actual Fermi energy", None),
                    defaults={"GRID": 1, "NE": 300, "ImE": 0.0})
    energy.copy_member("EMAX", is_optional=True)
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = InputParameters(
        "fscat",
        [CONTROL("FSCAT"), TAU, MODE, STRCONST, energy, TASK("FSCAT"), SITES],
        validators=_validate_fscat,
        info="Scattering amplitude matrices for three fixed backscattering geometries (ASA, serial).",
        description="All atomic types are evaluated. Initial directions are [001], [100], "
                    "and [111]; final directions are their opposites. No task-specific options are read. "
                    "The output energy axis is relative to EF in eV. FULLPOT, DMFT and LDA+U "
                    "are unsupported. SPR-KKR 9.7 has an erroneous spherical-harmonic index "
                    "and computes phase before guarding zero amplitude; numerical results are unverified."
    )
    schema._task_definition = TaskDefinition("FSCAT", executable="kkrgen", mpi=False, output_files=("{dataset}_GEO?_*.agr",))
    return schema
