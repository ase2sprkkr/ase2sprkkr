"""Paramagnetic current density from integrated Green functions, SPR-KKR 9.7."""
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters_definitions import InputParametersDefinition as InputParameters
from ._sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


def _validate_jorb(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    if parameters.ENERGY.GRID()[0] != "5":
        issues.append(DataValidityWarning("JORB expects an occupied-state integration contour; the E-INTEGRAL default is GRID=5."))
    if parameters.ENERGY.SPLITSS() or len(parameters.ENERGY.GRID()) > 1 or len(parameters.ENERGY.NE()) > 1:
        issues.append(DataValidityError("JORB forces NEPATH=1 and ORBCURDNS uses the first energy path; no split contour is supported."))
    if parameters.MODE.MODE() is not None:
        issues.append(DataValidityError("JORB uses relativistic wave functions and angular couplings; use FREL."))
    return issues or None


def input_parameters():
    energy = ENERGY(emax=(None, "Omitted upper bound delegates to EF", None),
                    defaults={"GRID": 5, "NE": 30, "ImE": 0.0})
    energy.copy_member("EMAX", is_optional=True)
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    return InputParameters(
        "jorb", [CONTROL("JORB"), TAU, MODE, STRCONST, energy, TASK("JORB"), SITES],
        executable="kkrgen", mpi=False, validators=_validate_jorb,
        info="Valence paramagnetic current density in the full-potential formalism.",
        description="No task-specific options are read. JORB becomes E-INTEGRAL internally. "
                    "The final ORBCURDNS call is skipped by the native MPI postprocessing guard. "
                    "Only the fully relativistic serial workflow is exposed. Numerical output is unverified.",
    )
