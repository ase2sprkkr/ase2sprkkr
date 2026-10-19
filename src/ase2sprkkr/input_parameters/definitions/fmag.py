"""Magnetic form factor from integrated densities, SPR-KKR 9.7."""
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters_definitions import InputParametersDefinition as InputParameters
from ..task_definition import TaskDefinition
from .sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


def _validate_fmag(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    if parameters.ENERGY.GRID()[0] != "5":
        issues.append(DataValidityWarning("FMAG expects an occupied-state integration contour; native FMAG check requests GRID=5."))
    if parameters.ENERGY.SPLITSS() or len(parameters.ENERGY.GRID()) > 1 or len(parameters.ENERGY.NE()) > 1:
        issues.append(DataValidityError("FMAG forces NEPATH=1 and DENSITY uses the first energy path; no split contour is supported."))
    if parameters.MODE.MODE() in ("NREL", "SREL"):
        issues.append(DataValidityWarning("Non-spin-polarized FMAG can divide by zero while normalizing spin/orbital form factors."))
    return issues or None


def input_parameters_definition():
    energy = ENERGY(emax=(None, "Omitted upper bound delegates to EF", None),
                    defaults={"GRID": 5, "NE": 30, "ImE": 0.0})
    energy.copy_member("EMAX", is_optional=True)
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = InputParameters(
        "fmag",
        [CONTROL("FMAG"), TAU, MODE, STRCONST, energy, TASK("FMAG"), SITES],
        validators=_validate_fmag,
        info="Magnetic form factor from charge, spin and orbital densities (ASA).",
        description="No task-specific form-factor options are read. FMAG becomes E-INTEGRAL internally. "
                    "The momentum mesh is fixed to 100 points. FULLPOT stops in GEN. "
                    "Normalization of a zero magnetic moment is unguarded in SPR-KKR. Numerical output is unverified."
    )
    schema._task_definition = TaskDefinition("FMAG", executable="kkrgen", mpi=True, output_file_templates=("{dataset}_fmag_*.agr", "{dataset}_dens_*.agr"))
    return schema
