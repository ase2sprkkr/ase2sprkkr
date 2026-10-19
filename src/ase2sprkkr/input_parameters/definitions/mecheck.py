"""Matrix-element diagnostics implemented by ME_TEST_11 in SPR-KKR 9.7."""

from ...common.grammar_types import Integer, Keyword, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters_definitions import InputParametersDefinition as InputParameters, InputValueDefinition as V
from ..task_definition import TaskDefinition
from ._file_manifests import resolve_file
from .sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


class MecheckInputParametersDefinition(InputParameters):
    """MECHECK definition with the files actually enabled in public 9.7."""

    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        if parameters.TASK.INITIALSTATE() == "CORE" and parameters.ENERGY.NE()[0] > 1:
            # K_CALC_ME=2 disables DATASET_check_*.dat in public 9.7.
            # Only non-zero matrix elements get one of these energy scans.
            files.extend(resolve_file(directory, f"MPZAZB_{polarization}_*-*-*.dat")
                         for polarization in ("-", "z", "+"))
        return tuple(files)


def _print_default(option):
    return option._container._get_root_container().CONTROL.PRINT()


def _validate_mecheck(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    task = parameters.TASK
    if parameters.MODE.MODE() is not None:
        issues.append(DataValidityError("MECHECK requires fully relativistic IREL=3."))
    if parameters.MODE.OP() in ("DMFT", "LDA+U"):
        issues.append(DataValidityError("MECHECK does not support DMFT or LDA+U (GEN)."))
    if task.NMEFORM.is_set():
        issues.append(DataValidityWarning(
            "SPR-KKR 9.7 overwrites TASK.NMEFORM with 2: only NAB and ADA are compared."))
    if task.INITIALSTATE() == "CORE":
        shell = task.CL()
        angular_momentum = "SPDF".index(shell[1])
        if angular_momentum >= int(shell[0]):
            issues.append(DataValidityError("TASK.CL must describe a shell with l < n."))
        nl = parameters.SITES.NL()
        if nl is not None and len(set(nl)) == 1 and angular_momentum >= int(nl[0]):
            issues.append(DataValidityError("The selected core orbital requires SITES.NL > l."))
        for name in ("ERYDI", "ERYDF"):
            if task[name].is_set():
                issues.append(DataValidityWarning(
                    f"TASK.{name} is ignored for INITIALSTATE=CORE; use CL and ENERGY."))
    else:
        if task.CL.is_set():
            issues.append(DataValidityWarning("TASK.CL is ignored for INITIALSTATE=BAND."))
        for name in ("ERYDI", "ERYDF"):
            if abs(complex(*task[name]())) < 1e-6:
                issues.append(DataValidityWarning(
                    f"SPR-KKR replaces TASK.{name} with its BAND default when its magnitude is <1e-6."))
        if any(parameters.ENERGY[name].is_set() for name in
               ("GRID", "NE", "ImE", "EMIN", "EMAX", "EMINEV", "EMAXEV")):
            issues.append(DataValidityWarning(
                "MECHECK BAND uses TASK.ERYDI and TASK.ERYDF; ENERGY does not select these energies."))
    return issues or None


def input_parameters_definition():
    energy = ENERGY(emax=(None, "Upper CORE final-state energy in Ry; omitted means the Fermi energy", None),
                    defaults={"GRID": 3, "NE": 1, "ImE": 0.01})
    energy.copy_member("EMAX", is_optional=True)
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = MecheckInputParametersDefinition(
        "mecheck",
        [CONTROL("MECHECK"), TAU, MODE, STRCONST, energy,
         TASK("MECHECK", [
             V("ITSEL", Integer(min=1), 1, is_always_added=False,
               info="One-based atomic-type index; native default 1."),
             V("PRINT", int, _print_default, is_required=False, is_always_added=False,
               info="Override CONTROL.PRINT in ME_TEST_11 (key PRINT, not IPRINT)."),
             V("INITIALSTATE", Keyword("BAND", "CORE"), "BAND", is_always_added=False),
             V("CL", Keyword(*(f"{n}{s}" for n in range(1, 7) for s in "SPDF")),
               "2P", is_always_added=False, info="Initial CORE shell; native default 2P."),
             V("ERYDI", SetOf(float, length=2), [0.5, 0.01], is_always_added=False,
               info="Initial BAND energy {real, imaginary} in Ry. Omitted or near-zero means {0.5,0.01}."),
             V("ERYDF", SetOf(float, length=2), [0.8, 0.01], is_always_added=False,
               info="Final BAND energy {real, imaginary} in Ry. Omitted or near-zero means {0.8,0.01}."),
             V("NONMAG", False, info="Suppress the magnetic field of the selected atomic type."),
             V("NMEFORM", Integer(min=1, max=3), 3, is_always_added=False,
               info="Native reader default 3, but SPR-KKR 9.7 unconditionally overwrites this with 2."),
         ]), SITES],
        validators=_validate_mecheck,
        info="Compare NAB and ADA matrix elements (fully relativistic, serial).",
        description="BAND uses two complex TASK energies; CORE uses CL and the ENERGY mesh. "
                    "Full-potential branches exist but are numerically unverified. "
                    "The compiled K_CALC_ME=2 disables the J-to-Z consistency test."
    )
    schema._task_definition = TaskDefinition("MECHECK", executable="kkrgen", mpi=False)
    return schema
