"""SOCPAR: single-site SOC and exchange parameters (SPR-KKR 9.7)."""
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters_definitions import InputParametersDefinition as InputParameters
from ..task_definition import TaskDefinition
from .sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


def _validate_socpar(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    if parameters.MODE.OP() in ("DMFT", "LDA+U"):
        issues.append(DataValidityError("SOCPAR does not support DMFT or LDA+U (GEN)."))
    if parameters.MODE.MODE() is not None:
        issues.append(DataValidityWarning(
            "SOCPAR uses relativistic wave-function channels and forces BS-SOC; "
            "non-fully-relativistic MODE settings are unverified."))
    if parameters.MODE.SOC.is_set():
        issues.append(DataValidityWarning("SOCPAR overwrites MODE.SOC scaling with SOCTL=1e-3."))
    if parameters.ENERGY.ImE() != 0:
        issues.append(DataValidityWarning("SOCPAR resets ENERGY.ImE to zero and uses real energies."))
    if parameters.ENERGY.SPLITSS():
        issues.append(DataValidityWarning("SOCPAR disables ENERGY.SPLITSS and uses only the first path."))
    if len(parameters.ENERGY.NE()) > 1 or len(parameters.ENERGY.GRID()) > 1:
        issues.append(DataValidityWarning("SOCPAR uses only the first ENERGY path."))
    if parameters.ENERGY.GRID()[0] != "1":
        issues.append(DataValidityWarning(
            "SOCPAR uses the real part of ETAB; GRID=1 supplies a direct real-energy scan. "
            "Complex contours are projected onto the real axis."))
    nl = parameters.SITES.NL()
    if nl is not None:
        if max(nl) < 2:
            issues.append(DataValidityError("SOCPAR needs SITES.NL>=2 to produce at least a p curve."))
        if max(nl) > 51:
            issues.append(DataValidityError("SOCPAR has room for at most 50 orbital legends (NL<=51)."))
    return issues or None


def input_parameters_definition():
    energy = ENERGY(emax=(None, "Upper energy in Ry; omitted means the actual Fermi energy", None),
                    defaults={"GRID": 1, "NE": 300, "ImE": 0.0})
    energy.copy_member("EMAX", is_optional=True)
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = InputParameters(
        "socpar",
        [CONTROL("SOCPAR"), TAU, MODE, STRCONST, energy, TASK("SOCPAR"), SITES],
        validators=_validate_socpar,
        info="Calculate single-site spin-orbit and exchange parameters (ASA, serial).",
        description="All atomic types and l=1..NL-1 are evaluated; no task-specific options are read. "
                    "SPR-KKR forces BS-SOC, suppresses the magnetic field and sets SOC scaling to 1e-3. "
                    "FULLPOT, DMFT and LDA+U are unsupported. The exchange-output eV label "
                    "does not match its unconverted Ry values in the 9.7 source."
    )
    schema._task_definition = TaskDefinition("SOCPAR", executable="kkrgen", mpi=False, output_file_templates=("{dataset}_soc_*.agr", "{dataset}_exc_*.agr"))
    return schema
