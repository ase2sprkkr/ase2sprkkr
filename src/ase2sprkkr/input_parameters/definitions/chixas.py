"""Field-induced XAS dichroism with externally supplied susceptibilities."""
import numpy as np
from ._filenames import dataset_adsi
from ...common.grammar_types import Integer, SetOf, real
from ...common.warnings import DataValidityError
from ._file_manifests import resolve_file, xray_shell_label
from ._sections import CONTROL, TAU, MODE, CPA, STRCONST, ENERGY, TASK, SITES
from ..input_parameters_definitions import InputParametersDefinition as InputParameters, InputValueDefinition as V, InputSectionDefinition as Section
from ..task_definition import TaskDefinition


class ChixasInputParametersDefinition(InputParameters):
    """CHIXAS definition with the source-defined XAS output families."""

    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        dataset = dataset_adsi(parameters)
        shell = xray_shell_label(parameters.TASK.CL())
        prefix = f"{dataset}_*"
        spectrum = f"{prefix}_{shell}"

        files.append(resolve_file(directory, f"{spectrum}.rat"))
        files.extend(resolve_file(directory, f"{prefix}_{kind}.pol")
                     for kind in ("ss", "ls", "sl", "ll"))
        files.extend(resolve_file(directory, f"{spectrum}_spec_{kind}")
                     for kind in ("0", "s", "o", "sr"))
        files.extend(resolve_file(directory, f"{spectrum}_{kind}.rat")
                     for kind in ("SSZ", "SSX", "SST", "LLZ", "LLX", "LLT", "TOT"))

        print_level = parameters.TASK.PRINT()
        if print_level is None:
            print_level = parameters.CONTROL.PRINT()
        if print_level >= 1:
            files.append(resolve_file(directory, f"{dataset}_*_GG.dat"))
        return tuple(files)


def _validate_chixas(p, _values, why):
    if why == "set":
        return None
    errors = []
    def error(message):
        errors.append(DataValidityError(message))
    if p.MODE.MODE() in ("NREL", "SREL", "SP-SREL"):
        error("This CHIXAS definition covers FREL only.")
    if p.SCF.FULLPOT():
        error("Native RSP does not support CHIXAS with FULLPOT.")
    if p.TAU.BZINT() != "POINTS":
        error("This CHIXAS definition requires BZINT=POINTS.")
    if p.ENERGY.GRID().tolist() != ["3"] or len(p.ENERGY.NE()) != 1 or p.ENERGY.NE()[0] < 2 or p.ENERGY.SPLITSS():
        error("Use one GRID=3 spectral mesh with NE>=2, without SPLITSS.")
    if not np.isfinite(p.ENERGY.ImE()) or p.ENERGY.ImE() <= 0:
        error("ImE must be finite and positive.")
    t = p.TASK
    spin, orbital, dos = t.CHISPN(), t.CHIVV(), t.NEF(all_values=True)
    if spin is None or orbital is None or dos is None:
        error("CHISPN, CHIVV and numbered NEF arrays are required; supply physical data for every type.")
    else:
        count = len(spin)
        if len(orbital) != count or len(dos) != count:
            error("CHISPN, CHIVV and NEF1..NEFn must cover the same contiguous type count.")
        if t.ITXRAY() > count:
            error("ITXRAY exceeds the supplied type count.")
        if not np.all(np.isfinite(spin)) or not np.all(np.isfinite(orbital)) or not np.all(np.isfinite(dos)):
            error("Susceptibilities and NEF must be finite")
        if np.any(dos < 0):
            error("NEF must be nonnegative.")

    return errors or None


def input_parameters_definition():
    energy = ENERGY(emax=(4.0, "Upper spectral energy in Ry; explicit ASE choice.", None),
                    defaults={"GRID": [3], "NE": [180], "ImE": 0.01})
    energy.copy_member("NE", type=SetOf(Integer(min=2), min_length=1, max_length=2))
    schema = ChixasInputParametersDefinition(
        "chixas",
        [CONTROL("CHIXAS"), TAU, MODE, CPA, STRCONST, energy,
            Section("SCF", [V("FULLPOT", False)], is_optional=True, is_expert=True),
            TASK("CHIXAS", add=[
                V("ITXRAY", Integer(min=1), 1, info="Absorbing atom type; native key is ITXRAY, not IT."),
                V("CL", str, "2P", info="Initial core shell."),
                V("NEF", SetOf(real, min_length=1), is_repeated="NUMBERED", is_optional=True,
                  info="NEF1, NEF2, ...: l-resolved DOS at EF for every potential type; must match native normalization."),
                V("CHISPN", SetOf(real, min_length=1), is_optional=True, info="Per-type spin susceptibility in 10^-6 cm^3/mol."),
                V("CHIVV", SetOf(real, min_length=1), is_optional=True, info="Per-type total orbital susceptibility in 10^-6 cm^3/mol; split equally between spins."),
                V("ABSUNIT", False, info="Use absolute absorption units."),
                V("PRINT", Integer(min=0), is_optional=True, info="Native CHIXAS print control."),
            ]), SITES],
        validators=_validate_chixas,
        info="Field-induced magnetic circular dichroism in XAS",
        description="Serial periodic ASA FREL CHIXAS with supplied NEF, CHISPN and CHIVV for every potential type. "
                    "Input-only validation cannot establish agreement with the potential or DOS normalization. "
                    "Spin-resolved orbital input, restart and TAU cache workflows are outside this audit."
    )
    schema._task_definition = TaskDefinition("CHIXAS", experimental=True, executable="kkrchi", mpi=False)
    return schema
