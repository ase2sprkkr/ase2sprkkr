"""2D plot windows for ASA and experimental FULLPOT PLY output."""
import numpy as np
from ...common.grammar_types import Integer, Keyword, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters_definitions import InputParametersDefinition as InputParameters, InputValueDefinition as V
from ..task_definition import TaskDefinition
from .sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


def _validate_plot2dpr(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    if parameters.MODE.OP() in ("DMFT", "LDA+U"):
        issues.append(DataValidityError("PLOT2DPR does not support DMFT or LDA+U (GEN)."))
    task = parameters.TASK
    a, b, d = task.AVEC(), task.BVEC(), task.DVEC()
    if not all(np.all(np.isfinite(v)) for v in (a, b, d)):
        issues.append(DataValidityError("PLOT2DPR window vectors must be finite."))
        return issues
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0 or np.linalg.norm(np.cross(a, b)) <= 1e-12 * na * nb:
        issues.append(DataValidityError("PLOT2DPR requires nonzero, nonparallel AVEC and BVEC."))
    else:
        if abs(np.dot(a, b)) > 1e-12:
            issues.append(DataValidityWarning(
                "SPR-KKR orthogonalizes BVEC against AVEC while preserving the length of BVEC."))
        if task.NB() is None and int(np.floor(task.NA() * nb / na + 0.5)) < 1:
            issues.append(DataValidityError("The native automatic NB rounds to zero; set NB>=1."))
    if task.TASK() == "PLOT2DRHO":
        issues.append(DataValidityWarning(
            "ASA PLOT2DRHO still writes potential, not density; FULLPOT uses the experimental PLY routine."))
    if parameters.ENERGY.NE()[0] > 1:
        issues.append(DataValidityWarning("PLOT2DPR plots stored radial data; ENERGY.NE is not the spatial resolution."))
    return issues or None


def input_parameters_definition():
    energy = ENERGY(emax=(None, "Unused upper energy; omitted retains the native Fermi-energy default", None),
                    defaults={"GRID": 1, "NE": 1, "ImE": 0.0})
    energy.copy_member("EMAX", is_optional=True)
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    task = TASK("PLOT2DPR", [
        V("AVEC", SetOf(float, length=3), [2.0, 0.0, 0.0], is_always_added=False,
          info="First plot-window edge in Cartesian coordinates, units ALAT."),
        V("BVEC", SetOf(float, length=3), [0.0, 2.0, 0.0], is_always_added=False,
          info="Second edge in units ALAT; orthogonalized against AVEC with its length preserved."),
        V("DVEC", SetOf(float, length=3), [0.0, 0.0, 0.0], is_always_added=False,
          info="Plot-window displacement from the common origin in units ALAT."),
        V("NA", Integer(min=1), 50, is_always_added=False,
          info="Number of intervals along AVEC; regular FULLPOT grid would have NA+1 points."),
        V("NB", Integer(min=1), is_optional=True,
          info="Intervals along BVEC; omitted means NINT(NA*|BVEC|/|AVEC|), not a fixed 50."),
    ])
    task.copy_member("TASK", type=Keyword("PLOT2DPR", "PLOT2DPOT", "PLOT2DRHO"), default_value="PLOT2DPR")
    schema = InputParameters(
        "plot2dpr",
        [CONTROL("PLOT2DPR"), TAU, MODE, STRCONST, energy, task, SITES],
        validators=_validate_plot2dpr,
        info="2D windows for stored potential/density data (serial).",
        description="Select PLOT2DPOT or PLOT2DRHO through TASK.TASK. The ASA routine writes "
                    "only potential PLY files and ignores the density selection. FULLPOT first calls "
                    "an experimental PLY routine which terminates the program before the older grid writer. "
                    "Window vectors are Cartesian in units ALAT, not fractional lattice coordinates. "
                    "DMFT and LDA+U are unsupported; numerical output is unverified."
    )
    schema._task_definition = TaskDefinition("PLOT2DPR", executable="kkrgen", mpi=False, output_file_templates=(
            "PLY-V-ASA-*.ply",
            "POT-ASA-LINES.ply",
            "PLY-V-FP-*.ply",
            "PLY-V-MT-*.ply",
            "PLY-FP-RWS-*.ply",
            "PLY-FP-RATM-*.ply",
            "PLY-V-IS.ply",
            "PLY-FP-CELL-*.ply",
            "PLY-FP-RMT.ply",
            "PLY-FP-CELL-BOUND-V.ply",
            "PLY-FP-CELL-BOUND-0.ply",
            "PLY-FP-RATM.ply",
            "PLY-FP-RWS.ply",
        ))
    return schema
