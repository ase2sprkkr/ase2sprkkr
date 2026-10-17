"""Native PLOT2DPR input selectors and their geometry."""
import numpy as np
from ._sections import CONTROL, TASK, TAU, ENERGY, SITES
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError
from ..input_parameters_definitions import InputParametersDefinition, InputSectionDefinition as Section, InputValueDefinition as V


_PLY_OUTPUTS = (
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
)


def _validate_plot2d(p, _container, why):
    if why == "set":
        return None
    errors = []
    a, b, d = p.TASK.AVEC(), p.TASK.BVEC(), p.TASK.DVEC()
    if not all(np.all(np.isfinite(v)) for v in (a, b, d)) or np.linalg.norm(a) == 0 or np.linalg.norm(b) == 0 or np.linalg.norm(np.cross(a, b)) <= 1e-12 * np.linalg.norm(a) * np.linalg.norm(b):
        errors.append(DataValidityError("Plot window needs finite, nonzero, nonparallel AVEC/BVEC and finite DVEC."))
    if p.TASK.TASK() == "PLOT2DRHO" and not p.SCF.FULLPOT():
        errors.append(DataValidityError("PLOT2DRHO needs FULLPOT; native ASA routine plots only potential."))
    return errors or None


def plot2d_parameters(task):
    t = TASK(task).copy(add=[
        V("AVEC", SetOf(float, length=3), [2., 0., 0.], info="Plot window vector in atomic units."),
        V("BVEC", SetOf(float, length=3), [0., 2., 0.], info="Window vector; native orthogonalizes it against AVEC, preserving length."),
        V("DVEC", SetOf(float, length=3), [0., 0., 0.], info="Plot window displacement in atomic units."),
        V("NA", Integer(min=1), 50, info="Subdivisions along AVEC."),
        V("NB", Integer(min=1), 50, info="Subdivisions along BVEC; explicitly written rather than native aspect-ratio fallback."),
    ])
    scf = Section("SCF", [V("FULLPOT", task == "PLOT2DRHO", is_always_added=task == "PLOT2DRHO")],
        is_optional=task != "PLOT2DRHO", is_expert=task != "PLOT2DRHO")
    return InputParametersDefinition(task.lower(), [CONTROL(task), TAU,
        ENERGY(emin=None, emax=None, defaults={"GRID": 0, "NE": 1, "ImE": 0.0}), t, SITES, scf],
        executable="kkrgen", mpi=False, validators=_validate_plot2d,
        output_file_templates=_PLY_OUTPUTS,
        info=f"{task}: native plane-window graphics from an existing potential. ASA supports potential only. "
             "FULLPOT invokes the native test PLY writer, produces both potential/density data regardless of the "
             "selector, then terminates before the ordinary plane-grid writer. This is not a complete selective plotting implementation.",
        description="Requires consistent stored data. FULLPOT=true does not convert an ASA potential. "
                    "The native test writer has fixed sampling details and output filenames; no numerical run certified here.")
