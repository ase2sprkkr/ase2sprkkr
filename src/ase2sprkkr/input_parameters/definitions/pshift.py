"""Phase-shift diagnostics, based on PSHIFT in SPR-KKR 9.7."""
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters_definitions import InputParametersDefinition as InputParameters
from ..task_definition import TaskDefinition, format_file_template
from ._file_manifests import resolve_file
from .sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


class PshiftInputParametersDefinition(InputParameters):
    """PSHIFT definition with branch-dependent native XMGR outputs."""

    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        dataset = format_file_template(parameters, '{dataset}_{adsi}')
        files.append(resolve_file(directory, f"{dataset}_pshift_*.agr"))
        if parameters.MODE.MODE() is None:
            files.append(resolve_file(directory, f"{dataset}_logdrv_*.agr"))
        else:
            files.append(resolve_file(directory, f"{dataset}_t-matrix_*.agr"))
        return tuple(files)


def _validate_pshift(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    if parameters.MODE.OP() in ("DMFT", "LDA+U"):
        issues.append(DataValidityError("PSHIFT does not support DMFT or LDA+U (GEN)."))
    energy = parameters.ENERGY
    if energy.ImE() != 0:
        issues.append(DataValidityWarning("PSHIFT resets ENERGY.ImE to zero and uses real energies."))
    if energy.SPLITSS():
        issues.append(DataValidityWarning("PSHIFT disables ENERGY.SPLITSS."))
    if len(energy.NE()) > 1 or len(energy.GRID()) > 1:
        issues.append(DataValidityWarning("PSHIFT uses only the first ENERGY path."))
    if energy.GRID()[0] != "1":
        issues.append(DataValidityWarning(
            "PSHIFT unwraps phases along the ordered real parts of ETAB. "
            "GRID=1 gives a direct ascending real-energy scan; other grids may repeat or reverse energies."))
    if parameters.MODE.MODE() is None and parameters.MODE.SOC.is_set():
        issues.append(DataValidityWarning(
            "Relativistic PSHIFT calls DIRBS directly without SOCTL; MODE.SOC does not scale its SOC."))
    nl = parameters.SITES.NL()
    if nl is not None:
        if len(set(nl)) > 1 and parameters.MODE.MODE() is not None:
            issues.append(DataValidityWarning(
                "Nonrelativistic PSHIFT fills only each type's local NL but writes global NL curves; "
                "nonuniform SITES.NL may expose uninitialized output entries."))
        if max(nl) > 5 and parameters.MODE.MODE() is None:
            issues.append(DataValidityWarning(
                "PSHIFT formats 2*j with I1; NL>5 produces invalid legends for j>=11/2."))
    return issues or None


def input_parameters_definition():
    energy = ENERGY(emax=(1.0, "Upper absolute energy in Ry; 1 Ry includes the resonance search window", None),
                    defaults={"GRID": 1, "NE": 300, "ImE": 0.0})
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = PshiftInputParametersDefinition(
        "pshift",
        [CONTROL("PSHIFT"), TAU, MODE, STRCONST, energy, TASK("PSHIFT"), SITES],
        validators=_validate_pshift,
        info="Calculate single-site phase shifts and resonance diagnostics (ASA, serial).",
        description="All atomic types are evaluated; no task-specific options are read. "
                    "Supports fully relativistic and NREL/SREL/SP-SREL branches. "
                    "Energies on the output axis are absolute eV, not relative to the Fermi energy. "
                    "FULLPOT, DMFT and LDA+U are unsupported. "
                    "Relativistic output includes logarithmic derivatives; nonrelativistic output includes t matrices."
    )
    schema._task_definition = TaskDefinition("PSHIFT", executable="kkrgen", mpi=False)
    return schema
