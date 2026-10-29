"""Dipole matrix-element plots, based on ME_PLOT in SPR-KKR 9.7.

BAND scans initial electron energies at fixed photon energy, or photon energies
at one initial electron energy. CORE scans final electron energies. The source
fixes m_j=1/2 and polarization q=0; neither is an input parameter.
"""

from ._filenames import dataset_adsi
from ...common.grammar_types import Integer, Keyword, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters_definitions import InputParametersDefinition as InputParameters, InputValueDefinition as V
from ..task_definition import TaskDefinition
from ._file_manifests import orbital_label, resolve_file
from ._sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


class MeplotInputParametersDefinition(InputParameters):
    """MEPLOT definition with its native state/self-energy XMGR outputs."""

    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        dataset = dataset_adsi(parameters)
        task = parameters.TASK

        if task.INITIALSTATE() == "CORE":
            state_label = task.CL()
        else:
            angular_momentum = task.L()
            state_label = "VB" if angular_momentum is None else orbital_label(angular_momentum)

        sefin = task.SEFIN()
        if sefin:
            state_label = f"{state_label}_SE{sefin}"
        files.append(resolve_file(directory, f"{dataset}_{state_label}_*.agr"))

        if sefin:
            files.extend((
                resolve_file(directory, f"{dataset}_SEfin_E_SE{sefin}_*.agr"),
                resolve_file(directory, f"{dataset}_SEfin_r_SE{sefin}_*.agr"),
            ))
        return tuple(files)


def _iprint_default(option):
    return option._container._get_root_container().CONTROL.PRINT()


def _validate_meplot(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    task = parameters.TASK
    core = task.INITIALSTATE() == "CORE"
    if parameters.MODE.MODE() is not None:
        issues.append(DataValidityWarning(
            "MEPLOT constructs relativistic kappa/j/m_j channels. "
            "Other MODE.MODE settings are not verified; use fully relativistic mode."))

    if core:
        shell = task.CL()
        angular_momentum = "SPDF".index(shell[1])
        if angular_momentum >= int(shell[0]):
            issues.append(DataValidityError("TASK.CL must describe a shell with l < n."))
        ignored = [name for name in ("L", "EPHOT", "NEPHOT", "EPHOTMIN", "EPHOTMAX")
                   if task[name].is_set()]
        if ignored:
            issues.append(DataValidityWarning(
                f"MEPLOT ignores TASK.{', TASK.'.join(ignored)} for INITIALSTATE=CORE."))
    else:
        angular_momentum = task.L()
        if task.CL.is_set():
            issues.append(DataValidityWarning("TASK.CL is ignored for INITIALSTATE=BAND."))
        nphot = task.NEPHOT()
        low, high = task.EPHOTMIN(), task.EPHOTMAX()
        if nphot > 1:
            if low is None or high is None:
                issues.append(DataValidityError(
                    "NEPHOT>1 requires both TASK.EPHOTMIN and TASK.EPHOTMAX."))
            elif high < low:
                issues.append(DataValidityError("TASK.EPHOTMAX must be >= TASK.EPHOTMIN."))
            if parameters.ENERGY.NE()[0] != 1:
                issues.append(DataValidityError(
                    "Photon scans require ENERGY.NE=1: ME_PLOT otherwise overwrites "
                    "each photon curve at every initial electron energy."))
            if task.EPHOT.is_set():
                issues.append(DataValidityWarning("TASK.EPHOT is ignored in a photon-energy scan."))
            if task.SEFIN() != 0:
                issues.append(DataValidityWarning("SPR-KKR disables TASK.SEFIN when NEPHOT>1."))
        elif low is not None or high is not None:
            issues.append(DataValidityError(
                "Use TASK.EPHOT for one photon energy; EPHOTMIN/EPHOTMAX require NEPHOT>1. "
                "The SPR-KKR range reader does not reliably handle partial ranges."))

    nl = parameters.SITES.NL()
    if nl is not None and len(set(nl)) == 1 and angular_momentum is not None:
        if angular_momentum >= int(nl[0]):
            issues.append(DataValidityError("The selected initial orbital requires SITES.NL > l."))
    return issues or None


def input_parameters_definition():
    # EMAX is left unset so SPR-KKR can use the actual Fermi energy.
    energy = ENERGY(emax=(None, "Upper electron energy in Ry; omitted means the Fermi energy", None),
                    defaults={"GRID": 3, "NE": 300, "ImE": 0.01})
    energy.copy_member("EMAX", is_optional=True)
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = MeplotInputParametersDefinition(
        "meplot",
        [CONTROL("MEPLOT"), TAU, MODE, STRCONST, energy,
         TASK("MEPLOT", [
             V("ITSEL", Integer(min=1), 1, is_always_added=False,
               info="One-based atomic-type index. SPR-KKR defaults to 1."),
             V("IPRINT", int, _iprint_default, is_required=False, is_always_added=False,
               info="Override CONTROL.PRINT inside ME_PLOT. The input key is IPRINT, not PRINT."),
             V("INITIALSTATE", Keyword("BAND", "CORE"), "BAND", is_always_added=False,
               info="BAND: ENERGY specifies initial states; CORE: ENERGY specifies final states."),
             V("L", Integer(min=0), is_optional=True,
               info="Select initial BAND orbital l. If omitted, all initial band orbitals are included."),
             V("CL", Keyword(*(f"{n}{shell}" for n in range(1, 7) for shell in "SPDF")),
               "2P", is_always_added=False, info="Initial CORE shell. Default: 2P."),
             V("NONMAG", False,
               info="Suppress the selected type's magnetic field and apply nonmagnetic dipole selection rules. "
                    "This TASK flag is distinct from CONTROL.NONMAG."),
             V("SEFIN", Integer(min=0, max=3), 0, is_always_added=False,
               info="Final-state self energy: 0 none; 1 Hedin-Lundqvist; 2 Dirac-Hara; "
                    "3 Dirac-Hara real and Hedin-Lundqvist imaginary part. Disabled in photon scans."),
             V("EPHOT", float, 100.0, is_always_added=False,
               info="Fixed photon energy in eV, for INITIALSTATE=BAND. Default: 100 eV."),
             # Writing even NEPHOT=1 selects the range branch in ME_PLOT.
             V("NEPHOT", Integer(min=1), 1, is_always_added=False,
               info="Number of photon energies. For >1 supply both bounds and ENERGY.NE=1."),
             V("EPHOTMIN", float, is_optional=True, info="Minimum photon energy in eV (NEPHOT>1)."),
             V("EPHOTMAX", float, is_optional=True, info="Maximum photon energy in eV (NEPHOT>1)."),
         ]), SITES],
        validators=_validate_meplot,
        info="Plot dipole matrix elements for core/band initial states (serial).",
        description="ME_PLOT fixes initial m_j=1/2 and polarization q=0. "
                    "Full-potential branches exist but have not been tested here. "
                    "Photon scans must use a single initial electron energy."
    )
    schema._task_definition = TaskDefinition("MEPLOT", experimental=True, executable="kkrgen", mpi=False)
    return schema
