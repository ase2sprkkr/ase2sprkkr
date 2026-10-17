"""Fully relativistic AES, audited against SPR-KKR 9.7."""
from ...common.grammar_types import Integer, Keyword
from ...common.warnings import DataValidityError
from ..input_parameters_definitions import InputValueDefinition as V
from ._file_manifests import (
    RequiredTauInputParametersDefinition,
    core_level_output_name,
    resolve_file,
)
from ._sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


class AESInputParametersDefinition(RequiredTauInputParametersDefinition):
    """AES definition with the native atom/core-shell spectrum filename."""

    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        name = core_level_output_name(
            parameters, "aes", type_separator="", shell_separator="."
        )
        files.append(resolve_file(directory, name))
        return tuple(files)


def _validate_aes(parameters, _values, why):
    if why == "set":
        return None
    issues = []
    if parameters.MODE.MODE() is not None:
        issues.append(DataValidityError("AES requires fully relativistic IREL=3."))
    if parameters.MODE.OP() in ("DMFT", "LDA+U"):
        issues.append(DataValidityError("AES does not support DMFT or LDA+U (GEN)."))
    shell = parameters.TASK.CL()
    if "SPDF".index(shell[1]) >= int(shell[0]):
        issues.append(DataValidityError("TASK.CL requires l < n."))
    if shell[1] == "F":
        issues.append(DataValidityError("AES f shells need 14 core states, but STRME has only 12 entries."))
    if max(parameters.ENERGY.NE()) < 3:
        issues.append(DataValidityError("AES needs capacity for at least three TAU energies for interpolation."))
    return issues or None


def input_parameters():
    energy = ENERGY(emax=(None, "Omitted bound delegates to EF", None),
                    defaults={"GRID": 3, "NE": 300, "ImE": .01})
    energy.copy_member("EMAX", is_optional=True)
    return AESInputParametersDefinition(
        "aes", [CONTROL("AES"), TAU, MODE, STRCONST, energy,
        TASK("AES", [
            V("IT", Integer(min=1), 1, is_always_added=False),
            V("CL", Keyword(*(f"{n}{s}" for n in range(1, 6) for s in "SPDF")), "2P", is_always_added=False),
            V("NEME", Integer(min=3), 10, is_always_added=False, info="Coarse matrix-element mesh; actual count is min(NEME, NE from TAU)."),
            V("PRINTME", False, info="Write matrix-element diagnostics."),
            V("PRINTAME", False, info="Print angular-matrix-element diagnostics."),
        ]), SITES], executable="kkrgen", mpi=False,
        cache_file_templates=("ME", *(f"ME{i}" for i in range(1, 13))),
        validators=_validate_aes,
        info="Fully relativistic Auger electron spectroscopy (ASA, NQ=1).",
        description="Requires an initial-state TAU directory on a uniform increasing mesh with at least three energies. "
                    "The source stops for NQ>1. ENERGY.NE must allocate enough capacity for the TAU file. "
                    "FULLPOT, DMFT and LDA+U are unsupported. Numerical output is unverified.",
    )
