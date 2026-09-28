"""BSF input-parameters definition."""

from functools import cache
from typing import Any, List, Optional

from ...common.generated_configuration_definitions import Length
from ...common.grammar_types import Keyword, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters import InputSection
from ...physics.lattice_data import bravais_number
from ..input_parameters_definitions import (
    InputParametersDefinition as InputParameters,
    InputValueDefinition as V,
)
from .sections import CONTROL, CPA, ENERGY, MODE, SITES, STRCONST, TASK, TAU


EK = "EK"
KK = "KK"

EK_TASK_ITEMS = ("NK", "KPATH", "NKDIR", "KE")
KK_TASK_ITEMS = ("NK1", "NK2", "K1", "K2")


def _mode_hint(parameters):
    name = getattr(parameters, "_requested_task_name", "").upper()
    if name == "BSFEK":
        return EK
    if name == "BSFKK":
        return KK
    return None


def bsf_mode(parameters):
    """Return the BSF mode selected by the number of energy points."""
    ne = parameters["ENERGY"]["NE"]()
    if ne is None or len(ne) == 0:
        return _mode_hint(parameters) or KK
    return EK if ne[0] > 1 else KK


def _root(container):
    return container._get_root_container()


def _ek_item(_definition, container):
    return bsf_mode(_root(container)) != KK


def _ek_vectors_item(_definition, container):
    return _ek_item(_definition, container) and container["KPATH"]() is None


def _kk_item(_definition, container):
    return bsf_mode(_root(container)) != EK


def _ka_item(definition, container):
    return _kk_item(definition, container) or _ek_vectors_item(definition, container)


def _ka_numbered(option):
    return bsf_mode(_root(option)) == EK


def _ne_default(option):
    return [200] if _mode_hint(_root(option)) == EK else [1]


def _emin_default(option):
    return -0.2 if bsf_mode(_root(option)) == EK else 0.7


def _emaxev_default(option):
    if bsf_mode(_root(option)) == KK:
        return option._container["EMINEV"]()
    return None


def _emax_default(option):
    if bsf_mode(_root(option)) == KK:
        return option._container["EMIN"]()
    return -1.0


def _validate_bsf(
    parameters: Any, _values: Any, why: str
) -> Optional[List[DataValidityWarning]]:
    issues = []
    mode = bsf_mode(parameters)
    hint = _mode_hint(parameters)
    if hint and hint != mode:
        issues.append(
            DataValidityWarning(
                f"ENERGY.NE selects BSF-{mode}, but BSF-{hint} was requested"
            )
        )

    task = parameters["TASK"]
    incompatible = KK_TASK_ITEMS if mode == EK else EK_TASK_ITEMS
    incompatible = [name for name in incompatible if task[name].is_set()]
    if incompatible:
        issues.append(
            DataValidityError(
                f"TASK.{', TASK.'.join(incompatible)} cannot be used in BSF-{mode} mode. "
                 "The mode is chosen by setting number of energies in ENERGY.NE."
            )
        )

    if mode == EK:
        if why != "set" and task["KPATH"]() is None and (
            task["KA"]() is None or task["KE"]() is None
        ):
            issues.append(
                DataValidityError(
                    "Please, specify either TASK.KPATH or TASK.KA and TASK.KE"
                )
            )
        nkdir = task["NKDIR"]()
        if nkdir is not None and nkdir > 9:
            issues.append(DataValidityError("TASK.NKDIR cannot be greater than 9"))
    elif why != "set":
        missing = [name for name in KK_TASK_ITEMS if task[name]() is None]
        if missing:
            issues.append(
                DataValidityError(
                    f"TASK.{', TASK.'.join(missing)} are required in BSFKK mode"
                )
            )
    return issues or None


kpaths = {
    0: "Default short path: Γ-X / Γ-H (Bravais 13, 14)",
    1: "Full standard high-symmetry path (4, 11, 12, 13, 14)",
    2: "Reduced standard high-symmetry path (4, 11, 12, 13, 14)",
    3: "Short standard high-symmetry path (4, 11, 12, 13, 14)",
    4: "Short path: Γ-Y-T-Z / Γ-M / Γ-X-M / Γ-X / Γ-H-N-Γ (4, 11, 12, 13, 14)",
    5: "Selected path: Γ-X / K-Γ / Γ-X,Y,Z / Γ-L / Γ-H (4, 11, 12, 13, 14)",
    6: "Selected directions: Γ-Y / Γ-X,Y,Z (4, 13, 14)",
    7: "Selected paths: X-Γ + Γ-Y (4)",
    10: "Alternative FCC path with L-K branch (13)",
}

kpaths_full = {
    4: {
        "name": "orthorhombic primitive",
        "paths": {
            1: ("full standard path",
                "Γ-Σ-X-G-U-A-Z-Λ-Γ-∆-Y-H-T-B-Z + X-D-S-C-Y + U-P-R-E-T + S-Q-T"),
            2: ("reduced standard path",
                "Γ-Σ-X-G-U-A-Z-Λ-Γ-∆-Y-H-T-B-Z"),
            3: ("short standard path",
                "Γ-Σ-X-G-U-A-Z-Λ-Γ"),
            4: ("secondary standard path",
                "Γ-∆-Y-H-T-B-Z"),
            5: ("Γ-X direction",
                "Γ-Σ-X"),
            6: ("Γ-Y direction",
                "Γ-∆-Y"),
            7: ("X-Γ and Γ-Y directions",
                "X-Σ-Γ + Γ-∆-Y"),
        },
    },

    11: {
        "name": "hexagonal primitive",
        "paths": {
            1: ("full standard path",
                "Γ-Σ-M-T'-K-T-Γ-∆-A-R-L-S'-H-S-A + M-U-L + K-P-H"),
            2: ("reduced standard path",
                "Γ-Σ-M-T'-K-T-Γ-∆-A-R-L-S'-H-S-A"),
            3: ("short standard path",
                "Γ-Σ-M-T'-K-T-Γ-∆-A"),
            4: ("Γ-M direction",
                "Γ-Σ-M"),
            5: ("K-Γ direction",
                "K-T-Γ"),
        },
    },

    12: {
        "name": "cubic primitive",
        "paths": {
            1: ("full standard path",
                "Γ-∆-X-Y-M-V-R-Λ-Γ-Σ-M"),
            2: ("reduced standard path",
                "Γ-∆-X-Y-M-V-R-Λ-Γ"),
            3: ("short standard path",
                "Γ-∆-X-Y-M-V-R"),
            4: ("Γ-X-M path",
                "Γ-∆-X-Y-M"),
            5: ("three principal Γ-axis directions",
                "Γ-∆-X + Γ-∆-Y + Γ-∆-Z"),
        },
    },

    13: {
        "name": "cubic face-centered",
        "paths": {
            0: ("default path",
                "Γ-∆-X"),
            1: ("full standard path",
                "X-∆-Γ-Λ-L-Q-W-N-K-Σ-Γ + L-M-U-S-X-Z-W-B-U"),
            2: ("reduced standard path",
                "X-∆-Γ-Λ-L-Q-W-N-K-Σ-Γ"),
            3: ("short standard path",
                "X-∆-Γ-Λ-L"),
            4: ("Γ-X direction",
                "Γ-∆-X"),
            5: ("Γ-L direction",
                "Γ-Λ-L"),
            6: ("three principal Γ-axis directions",
                "Γ-∆-X + Γ-∆-Y + Γ-∆-Z"),
            10: ("alternative full path with L-K branch",
                 "X-∆-Γ-Λ-L + L-?-K + W-N-K-Σ-Γ + L-M-U-S-X-Z-W-B-U"),
        },
    },

    14: {
        "name": "cubic body-centered",
        "paths": {
            0: ("default path",
                "Γ-D-H"),
            1: ("full standard path",
                "Γ-D-H-G-N-Σ-Γ-Λ-P-F-H + N-D-P"),
            2: ("reduced standard path",
                "Γ-D-H-G-N-Σ-Γ-Λ-P-F-H"),
            3: ("short standard path",
                "Γ-D-H-G-N-Σ-Γ-Λ-P"),
            4: ("Γ-H-N-Γ path",
                "Γ-D-H-G-N-Σ-Γ"),
            5: ("Γ-H direction",
                "Γ-D-H"),
            6: ("three principal Γ-axis directions",
                "Γ-∆-X + Γ-∆-Y + Γ-∆-Z"),
        },
    },
}

def kpaths_for_atoms(atoms_or_bravais):
    bravais = bravais_number(atoms_or_bravais)
    paths = kpaths_full.get(bravais, {}).get('paths', {})
    return {
        number: label for number, label in paths.items()
    }



def _task_definition():

    def kpaths_cli_help(kpaths_full=kpaths_full, width=88):
        lines = [ ' Valid values for the given Bravais lattice',
                  ' ------------------------------------------',
                  '  ']

        for bravais, data in kpaths_full.items():
            lines.append(f"{bravais:2}:  {data['name']}")

            for kpath, (description, path) in data["paths"].items():
                prefix = f"    {kpath:2}  {description}"
                path_indent = " " * 8

                if len(prefix) + 2 + len(path) <= width:
                    lines.append(f"{prefix:<40}  {path}")
                else:
                    lines.append(prefix)
                    lines.append(f"{path_indent}{path}")

            lines.append("")

        lines.append(" other Bravais lattices require manual k-path settings")
        return "\n".join(lines).rstrip()

    out = TASK(
        "BSF",
        add=[
            V("NK", 300, condition=_ek_item, info="total number of k-points"),
            V(
                "KPATH",
                Keyword(kpaths),
                condition=_ek_item,
                info="Predefined path in k-space",
                description=kpaths_cli_help(),
                is_optional=True,
            ),
            V(
                "NKDIR",
                Length(
                    "KA",
                    "KE",
                    default_values=[[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
                ),
                condition=_ek_vectors_item,
                info="Number of directions treated in k-spaces",
                is_optional=False,
                is_required=False,
            ),
            V(
                "KA",
                SetOf(float, length=3),
                condition=_ka_item,
                is_repeated=V.Repeated.NUMBERED_IF(_ka_numbered),
                info=(
                    "First k-vector segment in k-space in multiples of 2π/a and "
                    "rectangular coordinates with * = 1, ...,NKDIR"
                ),
                is_optional=False,
                is_required=False,
            ),
            V(
                "KE",
                SetOf(float, length=3),
                condition=_ek_vectors_item,
                is_repeated="NUMBERED",
                info=(
                    "Last k-vector segment in k-space in multiples of 2π/a and "
                    "rectangular coordinates with * = 1, ...,NKDIR"
                ),
                is_optional=False,
                is_required=False,
            ),
            V(
                "NK1",
                int,
                condition=_kk_item,
                info="number of k-points along k1",
                is_optional=True,
            ),
            V(
                "NK2",
                int,
                condition=_kk_item,
                info="number of k-points along k2",
                is_optional=True,
            ),
            V(
                "K1",
                SetOf(float, length=3),
                condition=_kk_item,
                is_optional=True,
                info="first k-vector to span a two-dimensional region in k-space",
            ),
            V(
                "K2",
                SetOf(float, length=3),
                condition=_kk_item,
                is_optional=True,
                info="second k-vector to span a two-dimensional region in k-space",
            ),
        ],
    )
    out['KPATH'].choices_for_atoms = kpaths_for_atoms

    return out


def _energy_definition():
    energy = ENERGY(
        emin=(-0.2, "The energy at which the BSF mesh starts", None),
        emax=(-1.0, "The energy at which the BSF mesh ends", None),
        defaults={"GRID": 3, "NE": _ne_default, "ImE": 0.001},
    )
    energy["EMIN"]._energy_default_value = _emin_default
    energy["EMAX"]._energy_default_value = _emax_default
    energy["EMAXEV"]._energy_default_value = _emaxev_default
    return energy


class BSFTaskSection(InputSection):
    def k_path_gui(self, atoms, *, parent=None):
        from ase2sprkkr.gui.k_path import k_path_gui

        out = k_path_gui(atoms, parent=parent)
        if out:
            self.KPATH.clear()
            self.set(out)


@cache
def bsf_definition():
    """Return the single definition shared by the BSFEK and BSFKK aliases."""
    out = InputParameters(
        "bsf",
        [
            CONTROL("BLOCHSF"),
            TAU,
            _task_definition(),
            _energy_definition(),
            MODE,
            STRCONST,
            CPA,
            SITES,
        ],
        executable="kkrgen",
        mpi=True,
        info="BSF - Bloch spectral functions in the E-K or K-K plane",
        validators=_validate_bsf,
    )
    out["TASK"].result_class = BSFTaskSection
    return out


input_parameters = bsf_definition
