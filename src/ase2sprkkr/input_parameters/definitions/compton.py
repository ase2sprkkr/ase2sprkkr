"""Standard split-contour FREL Compton momentum densities."""
import numpy as np
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from .sections import TASK, CONTROL, TAU, ENERGY, SITES, STRCONST, MODE, CPA
from ._file_manifests import RequiredTauInputParametersDefinition, resolve_file
from ..input_parameters_definitions import InputValueDefinition as V
from ..task_definition import TaskDefinition


class ComptonInputParametersDefinition(RequiredTauInputParametersDefinition):
    """COMPTON definition with native 1D/2D momentum-density outputs."""

    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        dataset = parameters.CONTROL.DATASET()
        task = parameters.TASK
        if task.MD2D() or task.EMD2D():
            suffix = "emd2d.dat" if task.EMD2D() else "mcp2d.dat"
            files.extend(resolve_file(directory, f"{dataset}_{part}_{suffix}")
                         for part in ("BND", "COR", "TOT"))
        else:
            files.extend(resolve_file(directory, f"{dataset}_{part}.mcp")
                         for part in ("BND", "COR", "TOT"))
            if parameters.CONTROL.PRINT() > 0:
                files.append(resolve_file(directory, f"{dataset}_data.mcp"))
        return tuple(files)


def _validate_compton(p, _values, why):
    if why == "set":
        return None
    issues = []
    def error(message):
        issues.append(DataValidityError(message))
    if p.MODE.MODE() in ("NREL", "SREL", "SP-SREL"):
        error("COMPTON requires FREL.")
    if p.ENERGY.GRID().tolist() != ["5", "5"] or len(p.ENERGY.NE()) != 2 or not p.ENERGY.SPLITSS() or not p.TASK.SPLITSS():
        error("This COMPTON workflow requires two GRID=5 contours and SPLITSS in ENERGY and TASK.")
    if p.TAU.BZINT() != "POINTS":
        error("This periodic COMPTON workflow requires BZINT=POINTS.")
    for name in ("PNMAX", "PPMAX"):
        value = p.TASK[name]()
        if not np.isfinite(value) or value < 0:
            error(f"{name} must be finite and nonnegative.")
    vectors = {name: p.TASK[name]() for name in ("PNVEC", "PPVEC1", "PPVEC2")}
    for name, value in vectors.items():
        if value is not None and (not np.all(np.isfinite(value)) or np.linalg.norm(value) == 0):
            error(f"{name} must be a finite nonzero vector.")
    normal, first, second = vectors.values()
    if normal is not None:
        if second is not None:
            error("PPVEC2 is ignored when PNVEC is present; omit PPVEC2.")
        pair = (normal, first)
    else:
        if (first is None) != (second is None):
            error("Without PNVEC, supply both PPVEC1 and PPVEC2 or neither.")
        pair = (first, second)
    if all(value is not None for value in pair):
        a, b = pair
        norms = np.linalg.norm(a) * np.linalg.norm(b)
        if norms and abs(np.dot(a, b)) > 1e-10 * norms:
            error("Momentum basis vectors must be perpendicular; native code does not orthogonalize them.")
    shifts = [p.TASK[name]() for name in ("SHFT", "SHFT_PP", "SHFT_PN")]
    if sum(shifts) > 1:
        error("Select at most one momentum shift flag.")
    if (shifts[0] or shifts[1]) and p.TASK.NPP() % 2 == 0:
        error("Shifted transverse grids require odd NPP to retain both endpoints.")
    if (shifts[0] or shifts[2]) and p.TASK.NPN() % 2 == 0:
        error("Shifted longitudinal grids require odd NPN to retain both endpoints.")
    if p.TASK.EMD2D() and p.TASK.PNMAX() != 0:
        error("Native EMD2D applies the longitudinal scaling twice; only PNMAX=0 is audited safely.")
    if p.TASK.MD2D() and p.TASK.EMD2D():
        error("Choose MD2D or EMD2D.")
    if p.TASK.PNMAX() == 0 and not p.TASK.EMD2D():
        issues.append(DataValidityWarning("PNMAX=0 collapses the longitudinal grid to one point; set PNMAX for a profile."))
    return issues or None


def input_parameters_definition():
    energy = ENERGY(defaults={"GRID": [5, 5], "NE": [30, 50], "ImE": 0.01})
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    energy.copy_member("SPLITSS", default_value=True, is_always_added=True)
    tau = TAU.copy(add=[V("TAUDIR", str, is_optional=True, info="Existing DATASET.tau.dir directory by default."),
                        V("TAUPATH", str, is_optional=True, info="Parent of TAUDIR; native default is current directory.")])
    tau.copy_member("BZINT", default_value="POINTS")
    members = [
        V("SPLITSS", True, is_always_added=True, info="Required here as well as in ENERGY: COMPTON rereads this flag."),
        V("MEMODE", Integer(min=1, max=2), 1, info="1 interpolates radial matrix elements; 2 evaluates each momentum."),
        V("NPMAX", Integer(min=2, max=10000), 250, info="Radial interpolation point count for MEMODE=1."),
        V("PNMAX", 0.0, info="Longitudinal momentum maximum in atomic units; zero gives one point."),
        V("PPMAX", 10.0, info="Transverse momentum maximum in atomic units."),
        V("NPN", Integer(min=1), 50), V("NPP", Integer(min=1), 50),
        V("USE_PEQUIV", Integer(min=0, max=2), 0, info="Momentum symmetry reuse mode."),
    ]
    members += [V(name, SetOf(float, min_length=3, max_length=3), is_optional=True)
                for name in ("PNVEC", "PPVEC1", "PPVEC2")]
    members += [V(name, False) for name in ("MD2D", "EMD2D", "SHFT", "SHFT_PP", "SHFT_PN", "INSIDE_SPHERE", "SQQ_RESOLVED")]
    schema = ComptonInputParametersDefinition(
        "compton",
        [CONTROL("COMPTON"), tau, MODE, STRCONST, CPA,
                           energy, TASK("COMPTON").copy(members), SITES],
        validators=_validate_compton,
        info="Compton momentum density and profiles",
        description="Periodic FREL split-contour workflow using an existing compatible TAU directory. "
                           "ONLY_CORE and USE_ZJOVZZ workflows are not audited. NONMAG is allowed. "
                           "Native code forces the ZJ wavefunction convention even for FULLPOT."
    )
    schema._task_definition = TaskDefinition("COMPTON", executable="kkrgen", mpi=True)
    return schema
