"""Radial potential, density and shape-function plots (SPR-KKR 9.7)."""
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters_definitions import InputParametersDefinition as InputParameters
from ..task_definition import TaskDefinition
from .sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


def _validate_plotprs(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    if parameters.MODE.OP() in ("DMFT", "LDA+U"):
        issues.append(DataValidityError("PLOTPRS does not support DMFT or LDA+U (GEN)."))
    if parameters.TASK.HFF():
        issues.append(DataValidityWarning("PLOTPRS plots stored radial data and does not run the hyperfine-field calculation loop."))
    if parameters.ENERGY.NE()[0] > 1:
        issues.append(DataValidityWarning(
            "PLOTPRS reads existing radial data; ENERGY.NE does not define a plot spectrum."))
    return issues or None


def input_parameters_definition():
    # INIT_MOD_ENERGY requires NEMAX>0 even though no energy loop is run.
    energy = ENERGY(emax=(None, "Unused upper energy; omitted retains the native Fermi-energy default", None),
                    defaults={"GRID": 1, "NE": 1, "ImE": 0.0})
    energy.copy_member("EMAX", is_optional=True)
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = InputParameters(
        "plotprs",
        [CONTROL("PLOTPRS"), TAU, MODE, STRCONST, energy, TASK("PLOTPRS"), SITES],
        validators=_validate_plotprs,
        info="Plot existing radial potentials, densities and full-potential shape functions (serial).",
        description="PLOTPRS enables all three plot groups. ASA writes spherical potentials and "
                    "densities; FULLPOT adds nonspherical components and shape functions per radial mesh. "
                    "The potential must already contain the radial data. This task does not create "
                    "a self-consistent density or new shape functions. No task-specific options are read. "
                    "DMFT and LDA+U are unsupported. Selective PLOTPOT/PLOTRHO/PLOTSFN task aliases "
                    "are not registered by this definition; standalone PLOTRHO has a native allocation bug."
    )
    schema._task_definition = TaskDefinition("PLOTPRS", executable="kkrgen", mpi=False, output_file_templates=(
            "{dataset}_V_*.agr",
            "{dataset}_VNS_*.agr",
            "{dataset}_RHO_*.agr",
            "{dataset}_RNS_*.agr",
            "{dataset}_SFN_M*.agr",
        ))
    return schema
