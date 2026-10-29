"""Prepare the s-like positron state used by POSANI."""
import numpy as np
from ...common.grammar_types import Keyword
from ...common.warnings import DataValidityError, DataValidityWarning
from ._sections import TASK, CONTROL, TAU, ENERGY, SITES, STRCONST, MODE, CPA
from ..input_parameters_definitions import InputParametersDefinition as InputParameters, InputValueDefinition as V, InputSectionDefinition as Section
from ..task_definition import TaskDefinition


def _validate_posaniprep(p, _values, why):
    if why == "set":
        return None
    issues = []
    if p.MODE.MODE() in ("NREL", "SREL", "SP-SREL"):
        issues.append(DataValidityError("POSANIPREP requires FREL."))
    if p.SCF.FULLPOT():
        issues.append(DataValidityError("This POSANIPREP definition audits the ASA workflow only."))
    if p.ENERGY.GRID().tolist() != ["0"] or p.ENERGY.NE().tolist() != [1] or p.ENERGY.SPLITSS():
        issues.append(DataValidityError("POSANIPREP requires one energy: GRID=0, NE=1, no SPLITSS."))
    if p.ENERGY.EMIN() is None or not np.isfinite(p.ENERGY.EMIN()) or p.ENERGY.EMINEV() is not None:
        issues.append(DataValidityError("Set a finite absolute EMIN for the chosen positron state, without EMINEV."))
    if p.TAU.BZINT() != "POINTS":
        issues.append(DataValidityError("This periodic POSANIPREP definition requires BZINT=POINTS."))
    issues.append(DataValidityWarning("Choose EMIN for the bottom of the positron band; the default -0.2 Ry does not determine it. "
                                     "Native preparation also requires 3D, KMROT=0 and NQ=NT for its TAUQQ export."))
    return issues


def input_parameters_definition():
    energy = ENERGY(defaults={"GRID": [0], "NE": [1], "ImE": 0.0})
    tau = TAU.copy(
        add=[
            TAU["BZINT"].copy(default_value="POINTS"),
        ],
    )
    schema = InputParameters(
        "posaniprep",
        [CONTROL("POSANIPREP"), tau, MODE, STRCONST, CPA, energy, Section("SCF", [V("FULLPOT", False)], is_optional=True, is_expert=True),
                       TASK("POSANIPREP", add=[
                           V("POSVC", Keyword({"BM": "Native BM correlation potential", "DR": "Native DR correlation potential", "IPM": "Independent-particle model"}), "BM",
                             info="ASA electron-positron correlation potential: BM, DR, or independent-particle IPM."),
                       ]), SITES],
        validators=_validate_posaniprep,
        info="Prepare positron wavefunction and scattering data for POSANI",
        description="Periodic 3D ASA FREL preparation at a user-chosen energy near the positron band bottom. "
                    "Uses an electronic potential and writes DATASET_pos_dat.pan. Native code keeps s_1/2 only. "
                    "Requires KMROT=0 and NQ=NT for the current export; an existing DATASET_positron.rho can suppress the energy loop and must not be reused as this task's cache."
    )
    schema._task_definition = TaskDefinition("POSANIPREP", experimental=True, executable="kkrgen", mpi=False, output_files=("{DATASET_ADSI}_pos_dat.pan",))
    return schema
