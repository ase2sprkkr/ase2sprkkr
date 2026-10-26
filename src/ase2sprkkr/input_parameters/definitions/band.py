"""KKRSPEC band calculation, distinct from BSF and EKREL."""
import numpy as np
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError
from ..task_definition import TaskDefinition
from ..tasks._spec_files import SPEC_INPUT_FILES, SPEC_OUTPUT_FILES
from .arpes import input_parameters_definition as arpes_parameters
from ._sections import CONTROL
from ..input_parameters_definitions import InputValueDefinition as V, InputSectionDefinition as Section


def _validate_band(p, _values, why):
    if why == "set" or p.TASK.INPVER() == 0:
        return None
    errors = []
    def error(message):
        errors.append(DataValidityError(message))
    if p.TASK.TASK() != "BAND":
        error("This dedicated definition requires TASK BAND.")
    if not np.isfinite(p.SPEC_BAND.EIG_IME()) or p.SPEC_BAND.EIG_IME() <= 0:
        error("EIG_IME must be finite and positive; native uses a strict inequality.")
    if p.ENERGY.SPLITSS() or len(p.ENERGY.NE()) != 1:
        error("Use one BAND energy mesh without SPLITSS.")
    for name in ("EWORK_EV", "IMV_INI_EV", "IMV_FIN_EV"):
        if not np.isfinite(p.ENERGY[name]()) or p.ENERGY[name]() < 0:
            error(f"{name} must be finite and nonnegative.")
    if any(p.SPEC_EL[name]() is not None for name in ("K3", "K4", "NK3", "NK4")):
        error("Native SPEC_INPUT does not read K3/K4/NK3/NK4.")
    if p.SPEC_EL.K1() is None:
        error("This BAND definition requires an explicit K1 momentum mesh.")
    for vector, count in (("K1", "NK1"), ("K2", "NK2")):
        values = p.SPEC_EL[vector]()
        if values is not None and (not np.all(np.isfinite(values)) or p.SPEC_EL[count]() is None):
            error(f"{vector} must be finite and requires {count}.")
    return errors or None


def input_parameters_definition():
    out = arpes_parameters().copy(
        new_schema=True,
        name="band", add=[CONTROL("BAND"), Section("SPEC_BAND", [
        V("NEIGMAX", Integer(min=1), 100000, info="Capacity for accepted band roots; native stops when exhausted."),
        V("EIG_IME", 0.01, info="Positive cutoff on abs(Im(kz)); explicitly supplied to prevent damping-dependent native fallback."),
        V("CALC_TRANS_PROB", False, info="Also calculate TRLEED transmission probabilities; advanced numerical workflow."),
    ])], validators=_validate_band, info="KKRSPEC Bloch band calculation",
        description="KKRSPEC BAND uses transfer-matrix Bloch eigenvalues and a kz imaginary-part filter. "
                    "It is not BSF or EKREL. Default momentum mesh is one k_parallel point; spectral range is relative to EF. "
                    "Transmission, many-body, full-potential and parallel branches are not numerically certified.")
    task = out["TASK"].copy(remove=["CRYS_VEC"], add=[V("CRYS_VECS", False, info="Presence disables primitive crystal vectors for Miller indices.")])
    task.copy_member("TASK", default_value="BAND")
    energy = out["ENERGY"].copy()
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    energy.copy_member("EWORK_EV", info="Work function in eV, used by the shared spectroscopic setup.")
    energy.copy_member("IMV_INI_EV", default_value=0.0)
    energy.copy_member("IMV_FIN_EV", default_value=0.0)
    electron = out["SPEC_EL"].copy()
    electron.copy_member("K1", default_value=[0., 0.], is_always_added=True)
    electron.copy_member("NK1", type=Integer(min=1), default_value=1, is_always_added=True)
    electron.copy_member("NK2", type=Integer(min=1))
    schema = out.copy(add=[task, energy, electron])
    schema.task_selectors = ("BAND",)
    schema._task_definition = TaskDefinition("BAND", experimental=True, executable="kkrspec", mpi=False, output_reader="SPEC",
        input_files=SPEC_INPUT_FILES, output_files=SPEC_OUTPUT_FILES,
    )
    return schema
