"""WFPLOT input parameters, checked against SPR-KKR 9.7 public sources.

WFPLOT plots radial core or band wave functions and runs only in serial mode.
The implementation in wfplot.f rejects full-potential calculations. BAND uses
only ETAB(1, 1), not a spectrum over the supplied energy mesh.
"""

from ...common.grammar_types import Integer, Keyword, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters_definitions import InputParametersDefinition as InputParameters, InputValueDefinition as V
from ..task_definition import TaskDefinition
from ._file_manifests import orbital_label, resolve_file
from ._sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


class WfplotInputParametersDefinition(InputParameters):
    """WFPLOT definition with native CORE/BAND wave-function plot outputs."""

    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        dataset = parameters.CONTROL.DATASET()
        task = parameters.TASK

        if task.STATE() == "CORE":
            files.append(resolve_file(directory, f"{dataset}_{task.CL()}_*.agr"))
            return tuple(files)

        angular_momentum = task.L()
        label = "?" if angular_momentum is None else orbital_label(angular_momentum)
        files.extend((
            resolve_file(directory, f"{dataset}_{label}_*.agr"),
            resolve_file(directory, f"{dataset}_WF_match_l{label}_*.agr"),
        ))
        return tuple(files)


def _print_default(option):
    """An omitted TASK.PRINT retains the current CONTROL.PRINT."""
    return option._container._get_root_container().CONTROL.PRINT()


def _validate_wfplot(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    task = parameters.TASK
    state = task.STATE()
    if parameters.MODE.OP() in ("DMFT", "LDA+U"):
        issues.append(DataValidityError("WFPLOT does not support MODE.OP=DMFT or LDA+U (GEN)."))
    if parameters.MODE.MODE() is not None:
        issues.append(DataValidityWarning("WFPLOT forces fully relativistic IREL=3; MODE.MODE is ignored."))

    nl = parameters.SITES.NL()
    uniform_nl = int(nl[0]) if nl is not None and len(set(nl)) == 1 else None
    if state == "CORE":
        core = task.CL()
        angular_momentum = "SPDF".index(core[1])
        if angular_momentum >= int(core[0]):
            issues.append(DataValidityError("TASK.CL must describe a shell with l < n."))
        if task.L.is_set():
            issues.append(DataValidityWarning("TASK.L is ignored for STATE=CORE; TASK.CL selects the shell."))
    else:
        angular_momentum = task.L()
        if angular_momentum is None and uniform_nl is not None:
            angular_momentum = uniform_nl - 1
        if task.CL.is_set():
            issues.append(DataValidityWarning("TASK.CL is ignored for STATE=BAND."))
        if uniform_nl is not None and angular_momentum is not None and angular_momentum >= uniform_nl:
            issues.append(DataValidityError("TASK.L must be smaller than SITES.NL."))
        if parameters.ENERGY.NE()[0] > 1:
            issues.append(DataValidityWarning("WFPLOT uses only the first energy point ETAB(1,1)."))

    mj = int(task.MJ().split("/")[0]) / 2.0
    if angular_momentum is not None and abs(mj) > angular_momentum + 0.5:
        issues.append(DataValidityWarning("SPR-KKR reduces TASK.MJ to sign(MJ)*(l+1/2)."))
    return issues or None


def input_parameters_definition():
    energy = ENERGY(emax="emin", defaults={"GRID": 3, "NE": 1, "ImE": 0.01})
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = WfplotInputParametersDefinition(
        "wfplot",
        [
            CONTROL("WFPLOT"),
            TAU,
            MODE,
            STRCONST,
            energy,
            TASK("WFPLOT", [
                # ITSEL=1 is an ASE2SPRKKR convenience. It MUST be written:
                # SECTION_SET_INTEGER('ITSEL', ITSEL, 9999, 1) requires it.
                V("ITSEL", Integer(min=1), 1, is_required=True,
                  info="One-based atomic-type index, not a lattice-site index. Required by WFPLOT."),
                V("STATE", Keyword("BAND", "CORE"), "BAND", is_always_added=False,
                  info="Plot band or core radial wave functions."),
                V("PRINT", int, _print_default, is_required=False, is_always_added=False,
                  info="Override CONTROL.PRINT inside WFPLOT."),
                V("L", Integer(min=0), is_optional=True,
                  info="BAND orbital angular momentum. If omitted, SPR-KKR uses NLQ(IQAT(1,ITSEL))-1."),
                V("CL", Keyword(*(f"{n}{shell}" for n in range(1, 7) for shell in "SPDF")),
                  "2P", is_always_added=False,
                  info="CORE shell, such as 1S or 2P. SPR-KKR defaults to 2P."),
                V("MJ", Keyword(*(f"{sign}{n}/2" for sign in ("", "+", "-") for n in range(1, 62, 2))),
                  "1/2", is_always_added=False,
                  info="Half-integer m_j, written as e.g. +3/2 or -1/2; decimals are not accepted."),
            ]),
            SITES,
        ],
        validators=_validate_wfplot,
        info="Plot radial core or band wave functions (ASA, fully relativistic, serial).",
        description="BAND plots only the first ENERGY point. CORE uses CL and MJ. "
                    "FULLPOT, DMFT and LDA+U are unsupported by SPR-KKR 9.7."
    )
    schema._task_definition = TaskDefinition("WFPLOT", experimental=True, executable="kkrgen", mpi=False)
    return schema
