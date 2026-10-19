"""Absorptive optical conductivity, periodic ASA FREL workflow."""
import numpy as np
from ...common.grammar_types import Integer
from ...common.warnings import DataValidityError
from .sections import CONTROL, TAU, MODE, CPA, STRCONST, ENERGY, TASK, SITES
from ..input_parameters_definitions import InputParametersDefinition as InputParameters, InputValueDefinition as V, InputSectionDefinition as Section
from ..task_definition import TaskDefinition


def _validate_optabs(p, _values, why):
    if why == "set":
        return None
    errors = []
    def error(message):
        errors.append(DataValidityError(message))
    if p.MODE.MODE() in ("NREL", "SREL", "SP-SREL"):
        error("This OPTABS definition covers FREL only.")
    if p.SCF.FULLPOT():
        error("Public RSP does not allow OPTABS with FULLPOT.")
    if p.TAU.BZINT() != "POINTS":
        error("This OPTABS definition requires BZINT=POINTS.")
    if p.ENERGY.GRID().tolist() != ["0"] or p.ENERGY.NE().tolist() != [1] or p.ENERGY.SPLITSS():
        error("Use GRID=0, NE=1 without SPLITSS for OPTABS initialization; TASK.NEB controls the actual grid.")
    t = p.TASK
    if not np.isfinite(t.TEMPEL()) or t.TEMPEL() <= 1:
        error("TEMPEL must be finite and greater than 1 K; native KBT is uninitialized otherwise.")
    if not np.isfinite(t.OMEGA2()) or t.OMEGA2() <= 0:
        error("OMEGA2 must be finite and positive in eV.")
    if t.OMEGA1() != 0:
        error("Use OMEGA1=0 for the native commensurate frequency mesh.")
    if t.NOMEGA() < t.NEB():
        error("NOMEGA must be at least NEB to prevent native OMEGATAB overflow.")
    if not np.isfinite(p.ENERGY.ImE()) or p.ENERGY.ImE() <= 0:
        error("ImE must be finite and positive.")
    if p.MODE.OP() not in (None, "NONE", "BROOKS"):
        error("Many-body optical workflows are outside this OPTABS definition.")
    return errors or None


def input_parameters_definition():
    schema = InputParameters(
        "optabs",
        [CONTROL("OPTABS"), TAU, MODE, CPA, STRCONST,
            ENERGY(emin=None, emax=None, defaults={"GRID": [0], "NE": [1], "ImE": 0.01}),
            Section("SCF", [V("FULLPOT", False)], is_optional=True, is_expert=True),
            TASK("OPTABS").copy([
                V("NEB", Integer(min=3), 30, info="Initial-state energy resolution; also actual frequency count before thermal extension."),
                V("NOMEGA", Integer(min=3), 100, info="Allocated frequency capacity; native later resets the actual count to NEB."),
                V("OMEGA1", 0.0, is_always_added=True, info="Use zero: native replaces it by OMEGA2/NEB; arbitrary offsets miss energy-pair matching."),
                V("OMEGA2", 8.0, info="Maximum photon energy in eV."),
                V("TEMPEL", 300.0, info="Electronic temperature in K; must exceed 1 K in this native version."),
                V("NKT", Integer(min=1), 6, info="Thermal energy extension in multiples of kT."),
                V("FRAMETET", float, is_optional=True, info="Response-frame Euler polar angle in degrees."),
                V("FRAMEPHI", float, is_optional=True, info="Response-frame Euler azimuth in degrees."),
                V("FRAMEGAM", float, is_optional=True, info="Response-frame third Euler angle in degrees."),
            ]), SITES],
        validators=_validate_optabs,
        info="Absorptive optical conductivity",
        description="Standard periodic ASA FREL OPTABS. Native task name becomes OPTICS-ABS. "
                    "TASK.NEB controls energy and photon grids; initialization ENERGY bounds are unused by OPT_ABS. "
                    "Native low-temperature KBT and frequency allocation require safeguards. MPI not certified here."
    )
    schema._task_definition = TaskDefinition("OPTABS", executable="kkrchi", mpi=False, output_file_templates=("sig_?????_ISPR??_NV.dat", "sig_?????_ISPR??_VC.dat"))
    return schema
