"""Native AIPES local/incoherent angle-integrated PES workflow."""
import numpy as np
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError
from ..task_definition import TaskDefinition
from ..tasks._spec_files import SPEC_INPUT_FILES, SPEC_OUTPUT_FILES
from .arpes import input_parameters_definition as arpes_parameters
from ._sections import CONTROL
from ._energy import energy_relative_default
from ..input_parameters_definitions import InputValueDefinition as V


def _validate_aipes(p, _values, why):
    if why == "set" or p.TASK.INPVER() == 0:
        return None
    errors = []
    def error(message):
        errors.append(DataValidityError(message))
    if p.TASK.TASK() != "AIPES":
        error("This dedicated definition requires TASK AIPES.")
    for name in ("THETA", "PHI"):
        values = p.SPEC_EL[name]()
        if not np.all(np.isfinite(values)) or values[0] != values[-1]:
            error("This AIPES definition uses fixed finite electron geometry, not an angular quadrature.")
    if p.SPEC_EL.NT() != 1 or p.SPEC_EL.NP() != 1:
        error("This AIPES definition requires NT=NP=1; these counts do not define an angular integral.")
    if any(p.SPEC_EL[name]() is not None for name in ("KA", "K1", "K2", "K3", "K4", "NK1", "NK2", "NK3", "NK4")):
        error("Momentum scans are outside this AIPES definition.")
    for name in ("EWORK_EV", "IMV_INI_EV", "IMV_FIN_EV"):
        value = p.ENERGY[name]()
        if not np.isfinite(value) or value < 0:
            error(f"{name} must be finite and nonnegative.")
    if not np.isfinite(p.SPEC_PH.EPHOT()) or p.SPEC_PH.EPHOT() <= 0:
        error("EPHOT must be finite and positive in eV.")
    if p.ENERGY.SPLITSS() or len(p.ENERGY.NE()) != 1:
        error("Use a single AIPES energy mesh without SPLITSS.")
    return errors or None


def input_parameters_definition():
    out = arpes_parameters().copy(
        new_schema = True,
        name="aipes", add=[CONTROL("AIPES")], validators=_validate_aipes,
        info="Native angle-integrated valence photoemission with surface contribution",
        description="AIPES is a native local/incoherent photoemission approximation plus surface contribution, "
                    "not numerical integration of an ARPES angular mesh. Fixed electron geometry is audited here; "
                    "converge layers, angular-momentum cutoff and optical potentials.")
    task = out["TASK"].copy(remove=["CRYS_VEC"], add=[V("CRYS_VECS", False, info="Native presence disables primitive crystal vectors for Miller indices.")])
    task.copy_member("TASK", default_value="AIPES")
    energy = out["ENERGY"].copy()
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    energy["EMAXEV"].default_value = energy_relative_default('EMAX', 0.0)
    energy.copy_member("EWORK_EV", info="Work function in eV.")
    for name in ("NT", "NP"):
        out["SPEC_EL"].copy_member(name, type=Integer(min=1))
    schema = out.copy(add=[task, energy])
    schema.task_selectors = ("AIPES",)
    schema._task_definition = TaskDefinition("AIPES", experimental=True, executable="kkrspec", mpi=False, output_reader="SPEC",
        input_files=SPEC_INPUT_FILES, output_files=SPEC_OUTPUT_FILES,
    )
    return schema
