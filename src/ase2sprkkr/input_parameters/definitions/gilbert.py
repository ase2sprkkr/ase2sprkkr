"""Gilbert task input parameters definition"""

import numpy as np
from ...common.grammar_types import real, Keyword, Integer
from ...common.warnings import DataValidityError
from ._file_manifests import TauCacheInputParametersDefinition, resolve_file
from ._sections import CONTROL, TAU, ENERGY, SITES, TASK, MODE, CPA, STRCONST
from ..input_parameters_definitions import InputValueDefinition as V
from ..task_definition import TaskDefinition


class GilbertInputParametersDefinition(TauCacheInputParametersDefinition):
    def input_files(self, parameters, directory):
        files = list(super().input_files(parameters, directory))
        if parameters.TASK.SETFLUCT() == "M_T" and parameters.TASK.FLUCTFIL():
            files.append(resolve_file(directory, parameters.TASK.FLUCTFIL()))
        return tuple(files)


def _validate_gilbert(p, _values, why):
    if why == "set":
        return None
    issues = []
    def error(message):
        issues.append(DataValidityError(message))
    if p.MODE.MODE() in ("NREL", "SREL", "SP-SREL") or p.CONTROL.NONMAG():
        error("This GILBERT definition requires magnetic FREL.")
    if p.ENERGY.GRID().tolist() != ["0"] or p.ENERGY.NE().tolist() != [1] or p.ENERGY.SPLITSS():
        error("GIL calculates at EF; use GRID=0, NE=1 without SPLITSS for initialization.")
    if p.TAU.BZINT() != "POINTS":
        error("This standard periodic GILBERT workflow requires BZINT=POINTS.")
    task = p.TASK
    mode = task.SETFLUCT()
    for name in ("TMPMIN", "TMPMAX", "DMAGDT", "TDEBYE"):
        value = task[name]()
        if value is not None and (not np.isfinite(value) or value < 0):
            error(f"{name} must be finite and nonnegative in this workflow.")
    if mode == "M_T":
        if not task.FLUCTFIL():
            error("SETFLUCT=M_T requires FLUCTFIL.")
        if any(task[name]() is not None for name in ("TMPMIN", "TMPMAX", "NTMP")):
            error("M_T reads its temperature mesh from FLUCTFIL; omit TMPMIN/TMPMAX/NTMP.")
    elif mode == "MLIN":
        if task.TMPMIN() is None:
            error("SETFLUCT=MLIN requires an explicit TMPMIN.")
        if task.NFPHI() is not None and task.NFPHI() > 1 and task.DMAGDT() is None:
            error("MLIN with NFPHI>1 requires DMAGDT.")
        if task.NTMP() is not None and task.TMPMAX() is None:
            error("MLIN NTMP without TMPMAX would be ignored.")
    elif any(task[name]() is not None for name in ("TMPMIN", "TMPMAX", "NTMP")):
        if not all(task[name]() is not None for name in ("TMPMIN", "TMPMAX", "NTMP")):
            error("Lattice-only temperatures require TMPMIN, TMPMAX and NTMP together.")
    if task.TMPMIN() is not None and task.TMPMAX() is not None and task.TMPMIN() > task.TMPMAX():
        error("TMPMAX must not be below TMPMIN.")
    if mode != "M_T" and task.FLUCTFIL() is not None:
        error("FLUCTFIL is only used by M_T in this definition.")
    if mode != "MLIN" and task.DMAGDT() is not None:
        error("DMAGDT is only used by MLIN.")
    return issues or None


def input_parameters_definition():
    schema = GilbertInputParametersDefinition(
        "gilbert",
        [
            CONTROL("GILBERT"),
            TAU,
            MODE,
            CPA,
            STRCONST,
            ENERGY(emin=None, emax=None, defaults={"GRID": 0, "NE": 1}),
            SITES.copy(defaults={"NL": 4}),
            TASK("GILBERT").copy(
                [
                    V("NTMP", Integer(min=1), is_required=False, info="Number of temperature points used for α(T)"),
                    V(
                        "SETFLUCT",
                        Keyword(
                            {
                                "MLIN": "Linear temperature mesh with optional linear M(T)=1-DMAGDT*T spin fluctuations; lattice vibrations can be combined.",
                                "M_T": "Temperature mesh and magnetization read from FLUCTFIL; spin fluctuations and lattice vibrations can be combined.",
                            }
                        ),
                        default_value=None,
                        is_required=False,
                        info="Finite temperature calculation mode",
                    ),
                    V(
                        "NVIBRA",
                        Integer(min=1),
                        is_required=False,
                        info="Number of directions for atomic displacements representing thermal lattice vibrations",
                    ),
                    V("TMPMIN", real, is_required=False, info="Lower limit of the temperature region (SETFLUCT=MLIN)."),
                    V("TMPMAX", real, is_required=False, info="Upper limit of the temperature region (SETFLUCT=MLIN)"),
                    V(
                        "FLUCTFIL",
                        str,
                        is_required=False,
                        info="Data ﬁle which contains the information about "
                        "temperature dependent magnetization, taken from the "
                        "experiment or Monte Carlo simulations",
                    ),
                    V("NFTET", Integer(min=1), is_required=False, info="Number of grid points specifying θ angle (SETFLUCT=M_T)"),
                    V("NFPHI", Integer(min=1), is_required=False, info="Azimuthal spin grid for M_T or MLIN."),
                    V("DMAGDT", real, is_optional=True, info="Slope of normalized M(T) in MLIN; required when NFPHI>1."),
                    V("TDEBYE", real, is_optional=True, info="Debye temperature in K for lattice displacements."),
                    V("INITSVEC", Integer(min=0, max=1), is_optional=True, info="0 uses native canonical displacement directions; 1 constructs directions."),
                ]
            ),
        ],
        validators=_validate_gilbert,
        description="Standard periodic magnetic FREL Gilbert damping at EF. Native GIL requires KMROT=0. "
                    "RDTAU is set by native INITVAR; compatible initialization TAU data may be required. "
                    "MCS, RDLM and NLCPA workflows are not audited here.",
        info="The Gilbert damping parameter calculation"
    )
    schema._task_definition = TaskDefinition("GILBERT", executable="kkrchi", output_files=(
            "{dataset}alfa_x.dat",
            "{dataset}alfa_T.dat",
            "{dataset}_alfa_T.agr",
        ))
    return schema


""" The Gilbert damping input parameters definition"""
