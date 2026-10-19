"""Dedicated SPLEED input, with incident energies in eV."""
import numpy as np
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError
from ..task_definition import TaskDefinition
from ..tasks._spec_files import SPEC_INPUT_FILES, SPEC_OUTPUT_FILES
from .arpes import input_parameters_definition as arpes_parameters
from .sections import CONTROL, ENERGY
from ..input_parameters_definitions import InputValueDefinition as V


def _validate_spleed(p, _values, why):
    if why == "set" or p.TASK.INPVER() == 0:
        return None
    errors = []
    def error(message):
        errors.append(DataValidityError(message))
    if p.TASK.TASK() != "SPLEED":
        error("This dedicated definition requires TASK SPLEED.")
    if p.ENERGY.EMIN.is_set() or p.ENERGY.EMAX.is_set():
        error("Use EMINEV/EMAXEV for incident energies; the SPLEED reader rereads these as eV, not relative to EF.")
    lo, hi = p.ENERGY.EMINEV(), p.ENERGY.EMAXEV()
    if lo is None or hi is None or not np.isfinite(lo) or not np.isfinite(hi) or lo <= 0 or hi < lo:
        error("Set finite positive incident EMINEV<=EMAXEV in eV.")
    if len(p.ENERGY.NE()) != 1 or p.ENERGY.NE()[0] < 2 or p.ENERGY.SPLITSS():
        error("Use a single SPLEED mesh with NE>=2 without SPLITSS; native one-point energy handling is inconsistent.")
    if any(p.SPEC_EL[name]() is not None for name in ("KA", "K1", "K2", "K3", "K4", "NK1", "NK2", "NK3", "NK4")):
        error("This SPLEED definition audits angular incidence only, not momentum scans.")
    for name in ("THETA", "PHI"):
        bounds = p.SPEC_EL[name]()
        if not np.all(np.isfinite(bounds)) or bounds[0] != bounds[-1]:
            error("This SPLEED definition uses fixed finite incidence angles.")
    if p.SPEC_EL.NT() != 1 or p.SPEC_EL.NP() != 1:
        error("Fixed-angle SPLEED requires NT=NP=1.")
    for name in ("EWORK_EV", "IMV_INI_EV", "IMV_FIN_EV"):
        value = p.ENERGY[name]()
        if not np.isfinite(value) or value < 0:
            error(f"{name} must be finite and nonnegative.")
    if p.SPEC_EL.POL0.is_set() or p.SPEC_EL.POL0L.is_set():
        error("This SPLEED definition uses POL_E; POL0/POL0L alternatives are not audited.")
    return errors or None


def input_parameters_definition():
    energy = ENERGY(emin=(None, "Incident energy lower bound", 11.0), emax=(None, "Incident energy upper bound", 60.0),
                    add=[V("EWORK_EV", 4.2, info="Work function in eV."), V("IMV_INI_EV", 0.2, info="Initial optical potential imaginary part in eV."),
                         V("IMV_FIN_EV", 0.2, info="Final optical potential imaginary part in eV.")], defaults={"GRID": 1, "NE": 200})
    energy.copy_member("NE", type=SetOf(Integer(min=2), min_length=1, max_length=2))
    energy.copy_member("EMINEV", info="Incident energy lower bound in eV; NOT relative to EF in SPLEED.")
    energy.copy_member("EMAXEV", info="Incident energy upper bound in eV; NOT relative to EF in SPLEED.")
    out = arpes_parameters().copy(name="spleed", remove=["SPEC_PH"], add=[CONTROL("SPLEED"), energy],
        validators=_validate_spleed, info="Spin-polarized low-energy electron diffraction",
        description="Dedicated fixed-angle SPLEED. Incident energies are read as eV. Native SPOL default 1 is retained. "
                    "Modern input requires NE>=2 here; momentum and alternative polarization workflows are not audited.")
    out["TASK"].copy_member("TASK", default_value="SPLEED")
    out["SPEC_EL"].copy_member("NT", type=Integer(min=1))
    out["SPEC_EL"].copy_member("NP", type=Integer(min=1))
    out["SPEC_STR"].copy_member("N_LAYER", default_value=1)
    out["SPEC_STR"].copy_member("NLAT_G_VEC", default_value=30)
    task = out["TASK"].copy(remove=["CRYS_VEC"], add=[V("CRYS_VECS", False, info="Native key: disables primitive crystal vectors for Miller indices when present.")])
    schema = out.copy(add=[task])
    schema.task_selectors = ("SPLEED",)
    schema._task_definition = TaskDefinition(
        "SPLEED", executable="kkrspec", mpi=False, output_reader="ARPES",
        input_file_templates=SPEC_INPUT_FILES, output_file_templates=SPEC_OUTPUT_FILES,
    )
    return schema
