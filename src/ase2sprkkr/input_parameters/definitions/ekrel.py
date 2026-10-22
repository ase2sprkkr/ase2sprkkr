"""Ordered-system E(k) and Fermi-surface search, SPR-KKR 9.7."""
import numpy as np
from ...common.grammar_types import Integer, Keyword, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters_definitions import InputParametersDefinition as InputParameters, InputValueDefinition as V
from ..task_definition import TaskDefinition
from .sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU
from .bsf import kpaths


def _validate_ekrel(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    task = parameters.TASK
    low, high = task.EMIN(), task.EMAX()
    if not np.isfinite(low) or not np.isfinite(high) or high <= low:
        issues.append(DataValidityError("EKREL requires finite TASK.EMIN < TASK.EMAX."))
    if task['F-SURF']():
        eslice = task.ESLICE()
        if eslice is not None and not low <= eslice <= high:
            issues.append(DataValidityError("EKREL requires EMIN <= ESLICE <= EMAX."))
        if task.NKDIR() or task.KPATH.is_set():
            issues.append(DataValidityWarning(
                "F-SURF samples the Brillouin-zone grid, not the selected k-path; KDIRTAB is still called first."))
    elif task.ESLICE.is_set():
        issues.append(DataValidityWarning("TASK.ESLICE is ignored without F-SURF."))
    count = task.NKDIR()
    ka = {i: task[f"KA{i}"]() for i in range(1, 10) if task[f"KA{i}"]() is not None}
    ke = {i: task[f"KE{i}"]() for i in range(1, 10) if task[f"KE{i}"]() is not None}
    if count:
        for index in range(1, count + 1):
            a, b = ka.get(index), ke.get(index)
            if a is None or b is None:
                issues.append(DataValidityError(f"NKDIR={count} requires both KA{index} and KE{index}."))
            elif not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)) or np.array_equal(a, b):
                issues.append(DataValidityError(f"KA{index}/KE{index} must be finite and define a nonzero segment."))
        if task.KPATH.is_set():
            issues.append(DataValidityWarning("Explicit NKDIR overrides TASK.KPATH."))
    if any(i < 1 or i > count for i in set(ka) | set(ke)):
        issues.append(DataValidityError("KA/KE indices must lie in 1..NKDIR; set NKDIR for a manual path."))
    if parameters.ENERGY.NE().tolist() != [1]:
        issues.append(DataValidityWarning("EKREL uses TASK.NE for its search; ENERGY.NE is only initialization."))
    return issues or None


def input_parameters_definition():
    energy = ENERGY(emax=(None, "Initialization-only upper energy; native fallback EF", None),
                    defaults={"GRID": 1, "NE": 1, "ImE": 0.0})
    energy.copy_member("EMAX", is_optional=True)
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = InputParameters(
        "ekrel",
        [CONTROL("EKREL"), TAU, MODE, STRCONST, energy,
        TASK("EKREL", [
            V("NK", Integer(min=2), 1000, is_always_added=False, info="Requested number of k-points; path/grid construction may change it."),
            V("EMIN", float, 0.01, is_always_added=False, info="Search lower energy in absolute Ry (TASK, not ENERGY)."),
            V("EMAX", float, 1.0, is_always_added=False, info="Search upper energy in absolute Ry."),
            V("NE", Integer(min=2), 501, is_always_added=False, info="Number of search energies; >=2 avoids division by NE-1."),
            V("F-SURF", False, info="Search Fermi-surface bands on a whole-Brillouin-zone grid."),
            V("ESLICE", float, is_optional=True, info="F-SURF slice in absolute Ry; omitted means the actual EF."),
            V("KPATH", Keyword(kpaths), 1, is_always_added=False,
              info="Predefined path; availability depends on Bravais lattice, as in BSF."),
            V("NKDIR", Integer(min=0, max=9), 0, is_always_added=False, info="Manual segments, each requires KAi and KEi; native maximum 9."),
        ] + [V(f"{name}{i}", SetOf(float, length=3), is_optional=True,
                info="Manual segment endpoint in Cartesian multiples of 2*pi/a.")
             for i in range(1, 10) for name in ("KA", "KE")]), SITES],
        validators=_validate_ekrel,
        info="E(k) bands or Fermi-surface search for ordered systems.",
        description="CPA alloys are rejected by SPR-KKR: use BSF for disorder. TASK sets the actual "
                    "energy search and k-path; ENERGY is used for initialization. KDIRTAB is called "
                    "even for F-SURF, so its lattice restrictions remain relevant. "
                    "MPI is implemented in the native k-loop; numerical output is unverified."
    )
    schema._task_definition = TaskDefinition("EKREL", executable="kkrgen", mpi=True, output_file_templates=("{dataset}_{adsi}.bnd",))
    return schema
