"""Shared native PLOTPRS selectors (not a separate task)."""
from ._sections import CONTROL, TASK, TAU, ENERGY, SITES
from ...common.warnings import DataValidityError
from ..input_parameters_definitions import InputParametersDefinition, InputSectionDefinition as Section, InputValueDefinition as V


def _validate_plot(p, _container, why):
    if why != "set" and p.TASK.TASK() == "PLOTSFN" and not p.SCF.FULLPOT():
        return [DataValidityError("PLOTSFN requires FULLPOT data; native ASA writes no shape functions.")]


def plot_parameters(task, info):
    scf = Section("SCF", [V("FULLPOT", task == "PLOTSFN", is_always_added=task == "PLOTSFN",
        info="PLOTSFN needs an existing FULLPOT potential with shape-function data; this flag does not generate it.")],
        is_optional=task != "PLOTSFN", is_expert=task != "PLOTSFN")
    return InputParametersDefinition(task.lower(), [CONTROL(task), TAU,
        ENERGY(emin=None, emax=None, defaults={"GRID": 0, "NE": 1, "ImE": 0.0}), TASK(task), SITES, scf],
        executable="kkrgen", mpi=False, validators=_validate_plot, info=info,
        description="Selective native PLOTPRS output from the existing potential, without an SCF or DOS calculation. "
                    "Native initialization maps this selector to internal TASK=PLOTPRS while retaining its output flag.")
