"""Integrate positron density and calculate the annihilation lifetime."""
import numpy as np
from ...common.grammar_types import Keyword, SetOf, Integer
from ...common.warnings import DataValidityError, DataValidityWarning
from .sections import TASK, CONTROL, TAU, ENERGY, SITES, STRCONST, MODE, CPA
from ..input_parameters_definitions import InputParametersDefinition as InputParameters, InputValueDefinition as V, InputSectionDefinition as Section
from ..task_definition import TaskDefinition


def _validate_posanitau(p, _values, why):
    if why == "set":
        return None
    issues = []
    if p.MODE.MODE() in ("NREL", "SREL", "SP-SREL"):
        issues.append(DataValidityError("POSANITAU requires FREL."))
    if p.SCF.FULLPOT():
        issues.append(DataValidityError("This POSANITAU definition audits the ASA workflow only."))
    if p.ENERGY.GRID().tolist() != ["5"] or len(p.ENERGY.NE()) != 1 or p.ENERGY.SPLITSS():
        issues.append(DataValidityError("POSANITAU requires one GRID=5 integration contour without SPLITSS."))
    if p.ENERGY.EMIN() is None or not np.isfinite(p.ENERGY.EMIN()) or p.ENERGY.EMINEV() is not None:
        issues.append(DataValidityError("Set a finite absolute EMIN for the positron integration, without EMINEV."))
    if p.ENERGY.EMAXEV() is not None or (p.ENERGY.EMAX() is not None and not np.isfinite(p.ENERGY.EMAX())):
        issues.append(DataValidityError("Use a finite absolute EMAX, without EMAXEV, for this workflow."))
    if p.TAU.BZINT() != "POINTS":
        issues.append(DataValidityError("This periodic POSANITAU definition requires BZINT=POINTS."))
    issues.append(DataValidityWarning("Choose integration bounds for the positron band; omitted EMAX uses the electronic EF, which need not be appropriate. "
                                     "Native POSANITAU requires 3D; POSVC does not change the fixed lifetime enhancement factor."))
    return issues


def input_parameters_definition():
    energy = ENERGY(emax=(None, "Upper integration bound in Ry; omitted native value is EF.", None),
                    defaults={"GRID": [5], "NE": [30], "ImE": 0.01})
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    tau = TAU.copy()
    tau.copy_member("BZINT", default_value="POINTS")
    schema = InputParameters(
        "posanitau",
        [CONTROL("POSANITAU"), tau, MODE, STRCONST, CPA, energy, Section("SCF", [V("FULLPOT", False)], is_optional=True, is_expert=True),
                       TASK("POSANITAU").copy([
                           V("POSVC", Keyword({"BM": "Native BM correlation potential", "DR": "Native DR correlation potential", "IPM": "Independent-particle model"}), "BM",
                             info="ASA electron-positron correlation potential: BM, DR, or independent-particle IPM."),
                       ]), SITES],
        validators=_validate_posanitau,
        info="Positron annihilation rate and lifetime",
        description="Periodic 3D ASA FREL lifetime from integrated positron density. "
                    "Uses an electronic potential and may reuse DATASET_positron.rho. "
                    "Integration bounds need physical verification. Lifetime enhancement is fixed independently of POSVC."
    )
    schema._task_definition = TaskDefinition("POSANITAU", experimental=True, executable="kkrgen", mpi=False, cache_files=("{dataset}_positron.rho",))
    return schema
