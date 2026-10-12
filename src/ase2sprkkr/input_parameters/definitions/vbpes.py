"""Angle-integrated valence photoemission in KKRGEN (SPR-KKR 9.7)."""
import numpy as np
from ...common.grammar_types import Integer, Keyword, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters_definitions import InputValueDefinition as V
from ._file_manifests import RequiredTauInputParametersDefinition
from ._sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


def _validate_vbpes(parameters, _values, why):
    if why == "set":
        return None
    task = parameters.TASK
    issues = []
    if parameters.MODE.MODE() is not None:
        issues.append(DataValidityError("VBPES requires fully relativistic IREL=3."))
    for name in ("QVEC", "QVECPH", "NSPIN"):
        value = task[name]()
        if value is not None and (not np.all(np.isfinite(value)) or np.linalg.norm(value) == 0):
            issues.append(DataValidityError(f"TASK.{name} must be finite and nonzero."))
    if task.GEOMETRY() in (3, 4):
        if task.TETPH() is not None:
            if task.PHIPH() is None:
                issues.append(DataValidityError("TETPH requires PHIPH."))
            if task.QVECPH() is not None:
                issues.append(DataValidityWarning("TETPH/PHIPH override QVECPH."))
        elif task.QVECPH() is None:
            issues.append(DataValidityError("Custom geometry requires TETPH/PHIPH or QVECPH."))
    elif task.QVEC() is not None and np.all(np.isfinite(task.QVEC())) and np.linalg.norm(task.QVEC()) != 0:
        issues.append(DataValidityWarning("Built-in geometry overwrites one QVEC component without renormalizing."))
    if task.USETFIN() and task.USETAUNN():
        issues.append(DataValidityWarning("USETFIN disables USETAUNN in VBPES."))
    if max(parameters.ENERGY.NE()) < 1:
        issues.append(DataValidityError("VBPES requires a nonempty energy mesh."))
    return issues or None


def input_parameters():
    energy = ENERGY(emax=(None, "Omitted upper bound uses EF", None),
                    defaults={"GRID": 3, "NE": 300, "ImE": 0.01},
                    add=[V("ERANGE", float, is_optional=True, info="Occupied energy window below EF, in eV; overrides EMIN.")])
    energy.copy_member("EMAX", is_optional=True)
    return RequiredTauInputParametersDefinition(
        "vbpes", [CONTROL("VBPES"), TAU, MODE, STRCONST, energy,
        TASK("VBPES", [
            V("EPHOT", float, 1253.6, is_always_added=False, info="Photon energy in eV; values rounding to 1/2 select 1253.6/1486.6 eV."),
            V("DEFEP0", float, 0.0, is_always_added=False, info="Read and printed in eV; not used to shift energies in VBPES."),
            V("ME", Keyword("ADA", "NAB"), "ADA", is_always_added=False),
            V("MODE", Integer(min=1, max=3), 2, is_always_added=False, info="Matrix-element mode; FULLPOT forces 2."),
            V("EWORK", float, is_optional=True, info="Work function in eV; only reported as EVAC in this routine."),
            V("GEOMETRY", Integer(min=0, max=4), 0, is_always_added=False),
            V("QVEC", SetOf(float, length=3), is_optional=True),
            V("QVECPH", SetOf(float, length=3), is_optional=True),
            V("TETPH", float, is_optional=True, info="Photon polar angle in degrees for custom geometry."),
            V("PHIPH", float, is_optional=True, info="Photon azimuth in degrees for custom geometry."),
            V("NSPIN", SetOf(float, length=3), is_optional=True, info="Spin direction; omitted means minus the resulting QHAT."),
            V("ImEFIN", float, is_optional=True, info="Final-state broadening in Ry; omitted means Im(ETAB(NE,1))."),
            V("USETAUNN", False, info="Use or calculate the final-state TAU directory."),
            V("USETFIN", False, info="Single-site final state; disables USETAUNN."),
        ]), SITES], executable="kkrgen", mpi=False,
        output_file_templates=("arpes.dat", "arpes_spol.dat", "spin.matrix.xps"),
        validators=_validate_vbpes,
        info="Angle-integrated valence-band photoemission (KKRGEN).",
        description="Requires the initial-state TAU directory. Full relativity is required. "
                    "This is a separate implementation from KKRSPEC ARPES. GRV stops as unimplemented. "
                    "Numerical output and FULLPOT have not been verified.",
    )
