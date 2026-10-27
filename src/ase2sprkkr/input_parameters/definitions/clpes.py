"""Core-level photoemission through native KKRSPEC (serial scope)."""
import numpy as np
from .arpes import input_parameters_definition as arpes_parameters
from ..task_definition import TaskDefinition
from ..tasks._spec_files import SPEC_INPUT_FILES, SPEC_OUTPUT_FILES
from ._sections import CONTROL
from ...common.grammar_types import Integer, Keyword, SetOf
from ...common.warnings import DataValidityError
from ..input_parameters_definitions import InputSectionDefinition as Section, InputValueDefinition as V


def _validate_clpes(p, _container, why):
    if why == "set":
        return None
    errors = []
    def error(message):
        errors.append(DataValidityError(message))
    if p.TASK.TASK() != "CLPES":
        error("This definition requires TASK CLPES.")
    if p.TASK.INPVER() == 0:
        return errors
    if p.SPEC_PH.NPOL() is not None and p.SPEC_PH.NPOL() not in ("1", "2"):
        error("CLPES does not implement NPOL=0 or 3.")
    if p.ENERGY.SPLITSS() or len(p.ENERGY.NE()) != 1:
        error("CLPES requires a single energy mesh without SPLITSS.")
    if p.ENERGY.NE()[0] == 1 and p.SPEC_CL.IE_ICST() > 4 * "SPDF".index(p.SPEC_CL.CL()[1]) + 2:
        error("IE_ICST exceeds the selected core shell state count.")
    if p.SPEC_EL.NT() != 1 or p.SPEC_EL.NP() != 1:
        error("This CLPES definition audits fixed electron angles only.")
    for name in ("THETA", "PHI"):
        v = p.SPEC_EL[name]()
        if not np.all(np.isfinite(v)) or v[0] != v[-1]:
            error("Use fixed finite electron angles for CLPES.")
    if any(p.SPEC_EL[n]() is not None for n in ("KA", "K1", "K2", "K3", "K4", "NK1", "NK2", "NK3", "NK4")):
        error("Momentum scans are outside this CLPES definition.")
    for name in ("EWORK_EV", "IMV_INI_EV", "IMV_FIN_EV"):
        if not np.isfinite(p.ENERGY[name]()) or p.ENERGY[name]() < 0:
            error(f"{name} must be finite and nonnegative.")
    if not np.isfinite(p.SPEC_PH.EPHOT()) or p.SPEC_PH.EPHOT() <= 0:
        error("EPHOT must be finite and positive; also choose it above the core emission threshold.")
    return errors or None


def input_parameters_definition():
    out = arpes_parameters()
    task = out["TASK"].copy(remove=["CRYS_VEC"], add=[
        out["TASK"]["TASK"].copy(type=Keyword("CLPES"), default_value="CLPES"),
        V("CRYS_VECS", False),
    ])
    energy = out["ENERGY"]
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    photon = out["SPEC_PH"]
    photon.copy_member("EPHOT", default_value=1000.0)
    for name in ("EMIN", "EMAX", "EMINEV", "EMAXEV"):
        energy.copy_member(name, info="Initialization bound; CLPES replaces the spectral window using core energies +/-3 eV.")
    cl = Section("SPEC_CL", [
        V("IT", Integer(min=1), 1, info="Potential type; native selector also includes sites of the same atomic number."),
        V("CL", Keyword(*[f"{n}{l}" for n in range(1, 8) for i, l in enumerate("SPDF") if i < n]), "2P", is_always_added=True),
        V("IE_ICST", Integer(min=1), 1, info="Core state index, used only for NE=1."),
    ])
    schema = out.copy(
        name="clpes",
        validators=_validate_clpes,
        info="Core-level photoemission and photoelectron diffraction",
        description="Ordered potential only: native CLPES stops for CPA. Energy window is generated from core levels +/-3 eV; NE=1 selects IE_ICST. Photon energy must exceed the emission threshold.",
        add=[CONTROL("CLPES"), task, energy, photon, cl],
    )
    schema.task_selectors = ("CLPES",)
    schema._task_definition = TaskDefinition("CLPES", experimental=True, executable="kkrspec", mpi=False, output_reader="ARPES",
        input_files=SPEC_INPUT_FILES,
        output_files=(*SPEC_OUTPUT_FILES, "{dataset}_corewf")
    )
    return schema
