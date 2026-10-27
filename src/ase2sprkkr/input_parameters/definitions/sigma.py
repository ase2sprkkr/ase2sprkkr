"""Periodic ASA FREL DC conductivity, surface or sea plus surface."""
import numpy as np
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError
from ._sections import CONTROL, TAU, MODE, CPA, STRCONST, ENERGY, TASK, SITES
from ..input_parameters_definitions import InputParametersDefinition as InputParameters, InputValueDefinition as V, InputSectionDefinition as Section
from ..task_definition import TaskDefinition


def _validate_sigma(p, _values, why):
    if why == "set":
        return None
    errors = []
    def error(message):
        errors.append(DataValidityError(message))
    if p.MODE.MODE() in ("NREL", "SREL", "SP-SREL"):
        error("This SIGMA definition covers FREL only.")
    if p.SCF.FULLPOT():
        error("This SIGMA definition audits the ASA workflow only.")
    if p.TAU.BZINT() != "POINTS":
        error("This periodic SIGMA definition requires BZINT=POINTS.")
    ne = p.ENERGY.NE().tolist()
    grid = p.ENERGY.GRID().tolist()
    if len(ne) != 1 or p.ENERGY.SPLITSS() or not ((ne == [1] and grid == ["0"]) or (len(ne) == 1 and ne[0] > 1 and grid == ["5"])):
        error("Use GRID=0, NE=1 for Fermi surface only, or one GRID=5 contour with NE>1 for sea plus surface; no SPLITSS.")
    if not np.isfinite(p.ENERGY.ImE()) or p.ENERGY.ImE() <= 0:
        error("SIGMA requires positive ImE in this workflow.")
    projectors = p.TASK.SIGSPROJ()
    if projectors is not None:
        values = projectors.tolist()
        if values[0] not in (1, 9):
            error("The first SIGSPROJ must be 1 or 9, as required by native SIG.")
        if len(values) != len(set(values)):
            error("SIGSPROJ must not contain duplicates.")
    if p.MODE.OP() not in (None, "NONE", "BROOKS"):
        error("Many-body SIGMA workflows are outside this definition.")
    return errors or None


def input_parameters_definition():
    energy = ENERGY(defaults={"GRID": [0], "NE": [1], "ImE": 0.01})
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = InputParameters(
        "sigma",
        [CONTROL("SIGMA"), TAU, MODE, CPA, STRCONST, energy,
            Section("SCF", [V("FULLPOT", False)], is_optional=True, is_expert=True),
            TASK("SIGMA", add=[
                V("SIGSPROJ", SetOf(Integer(min=1, max=16), min_length=1, max_length=16), is_optional=True,
                  info="Response operator indices; first must be 1 or 9. Omitted native FREL selection is 1..16."),
                V("SP_NO_ALPHA", False, info="Disable alpha part of the spin-polarization operator."),
                V("SP_NO_NABLA", False, info="Disable nabla part of the spin-polarization operator."),
                V("NOSYMSIG", False, info="Disable response and global native symmetry."),
                V("FRAMETET", float, is_optional=True, info="Response-frame Euler polar angle in degrees."),
                V("FRAMEPHI", float, is_optional=True, info="Response-frame Euler azimuth in degrees."),
                V("FRAMEGAM", float, is_optional=True, info="Response-frame third Euler angle in degrees."),
            ]), SITES],
        validators=_validate_sigma,
        info="DC conductivity and selected response tensors",
        description="Standard periodic ASA FREL Kubo-Streda SIGMA. Default NE=1 omits the Fermi-sea integral; "
                    "use GRID=5 and NE>1 to include it. Native surface point is EF+i*ImE. "
                    "Layer, NLCPA, PROJECT, thermal and many-body workflows are not audited here; MPI not certified."
    )
    schema._task_definition = TaskDefinition("SIGMA", experimental=True, executable="kkrchi", mpi=False, output_files=("sigma?????.dat",))
    return schema
