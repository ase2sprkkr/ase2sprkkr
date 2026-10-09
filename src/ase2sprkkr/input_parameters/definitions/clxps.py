"""Core-level XPS input, audited against CLXPS in SPR-KKR 9.7."""
import numpy as np
from ...common.grammar_types import Integer, Keyword, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters_definitions import InputParametersDefinition as InputParameters, InputValueDefinition as V
from ._sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


def _validate_clxps(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    task = parameters.TASK
    if parameters.MODE.MODE() is not None:
        issues.append(DataValidityError("CLXPS requires fully relativistic IREL=3."))
    if parameters.MODE.OP() in ("DMFT", "LDA+U"):
        issues.append(DataValidityError("CLXPS does not support DMFT or LDA+U (GEN)."))
    shell = task.CL()
    l = "SPDF".index(shell[1])
    if l >= int(shell[0]):
        issues.append(DataValidityError("TASK.CL must describe a shell with l < n."))
    if max(parameters.ENERGY.NE()) < 4*l+2:
        issues.append(DataValidityError("ENERGY.NE must allocate at least 4*l+2 core-state energies for CLXPS."))
    if abs(task.ImE()) < 1e-8:
        issues.append(DataValidityWarning("CLXPS replaces near-zero TASK.ImE with 0.01 Ry."))
    present = []
    for i in range(11):
        values = [task[f"{name}{i}"]() for name in ("ANGQ", "ANGE", "ANGK")]
        if any(v is not None for v in values):
            present.append(i)
            if any(v is None for v in values):
                issues.append(DataValidityError(f"Geometry {i} requires ANGQ{i}, ANGE{i} and ANGK{i}."))
            elif not all(np.all(np.isfinite(v)) for v in values):
                issues.append(DataValidityError("CLXPS geometry angles must be finite."))
    if 0 in present:
        if len(present) > 1:
            issues.append(DataValidityWarning("Geometry 0 overrides numbered geometries 1..10."))
    elif present and present != list(range(1, max(present)+1)):
        issues.append(DataValidityError("CLXPS geometries must be consecutive starting at 1."))
    if task.ANGRES() and not present:
        issues.append(DataValidityError("ANGRES requires at least one complete geometry."))
    if not task.ANGRES() and len(present) > 1:
        issues.append(DataValidityWarning("Angle-integrated CLXPS uses only one geometry."))
    return issues or None


def input_parameters():
    # Allocate enough space even for 4f (14 states); native CLXPS default is 18.
    energy = ENERGY(emax=(None, "Initialization bound; actual final energies come from core levels plus EPHOT", None),
                    defaults={"GRID": 3, "NE": 18, "ImE": 0.01})
    energy.copy_member("EMAX", is_optional=True)
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    values = [
        V("ME", Keyword("ADA"), "ADA", is_always_added=False, info="Only ADA reaches the implemented matrix-element branch."),
        V("IT", Integer(min=1), 1, is_always_added=False, info="One-based atomic type index."),
        V("CL", Keyword(*(f"{n}{s}" for n in range(1, 7) for s in "SPDF")), "2P", is_always_added=False),
        V("EPHOT", float, 1253.6, is_always_added=False, info="Photon energy in eV; final states must lie above EF."),
        V("ImE", float, 0.01, is_always_added=False, info="Final-state imaginary energy in Ry; near-zero values become 0.01."),
        V("MECHECK", False, info="Print diagnostic matrix elements."),
        V("LINPOL", False, info="Linear polarization; omitted means circular."),
        V("USETSS", False, info="Use the single-site t matrix instead of the scattering-path operator."),
        V("ANGRES", False, info="Angle-resolved case; supply a complete geometry."),
    ] + [V(f"{name}{i}", SetOf(float, length=2), is_optional=True,
           info="Polar and azimuthal angles in degrees; photon q, polarization e or electron k.")
         for i in range(11) for name in ("ANGQ", "ANGE", "ANGK")]
    return InputParameters(
        "clxps", [CONTROL("CLXPS"), TAU, MODE, STRCONST, energy, TASK("CLXPS", values), SITES],
        executable="kkrgen", mpi=False, validators=_validate_clxps,
        info="Core-level XPS (ASA, fully relativistic).",
        description="TASK sets core shell, photon energy and optional angular geometries. "
                    "Final electron energies are core eigenenergies + EPHOT, not an ENERGY scan. "
                    "Geometries use index 0 for one overriding geometry, or consecutive indices 1..10. "
                    "FULLPOT, DMFT and LDA+U are unsupported; numerical output is unverified.",
    )
