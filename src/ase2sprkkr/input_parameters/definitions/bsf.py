"""BSF input-parameters definition."""

from functools import cache
from typing import Any, Optional

from ...common.generated_configuration_definitions import Length
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters import InputSection
from ..input_parameters_definitions import (
    InputParametersDefinition as InputParameters,
    InputValueDefinition as V,
)
from .sections import CONTROL, ENERGY, MODE, SITES, STRCONST, TASK, TAU


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
    if option._container["EMINEV"]() is not None:
        return None
    return -0.2 if bsf_mode(_root(option)) == EK else 0.7


def _emaxev_default(option):
    if bsf_mode(_root(option)) == KK:
        return option._container["EMINEV"]()
    return None


def _emax_default(option):
    if option._container["EMAXEV"]() is not None:
        return None
    if bsf_mode(_root(option)) == KK:
        return option._container["EMIN"]()
    return -1.0


def _validate_bsf(
    parameters: Any, _values: Any, why: str
) -> Optional[list[DataValidityWarning]]:
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


def _validate_kpath(option, _container, _why):
    value = option()
    if value is not None and value not in (*range(1, 8), 10):
        return DataValidityError("TASK.KPATH has to be one of 1-7 or 10")


def _task_definition():
    return TASK(
        "BSF",
        add=[
            V("NK", 300, condition=_ek_item, info="total number of k-points"),
            V(
                "KPATH",
                Integer(min=1, max=10),
                condition=_ek_item,
                validators=_validate_kpath,
                info="Predefined path in k-space",
                description="""
Bravais-lattice KPATH path
==========================
orb 1  Γ-Σ-X-G-U-A-Z-Λ-Γ-∆-Y-H-T-B-Z
       + X-D-S-C-Y + U-P-R-E-T + S-Q-T
    2  Γ-Σ-X-G-U-A-Z-Λ-Γ-∆-Y-H-T-B-Z
    3  Γ-Σ-X-G-U-A-Z-Λ-Γ
    4  Γ-∆-Y-H-T-B-Z
hex 1  Γ-Σ-M-T’-K-T-Γ-∆-A-R-L-S’-H-S-A
       + M-U-L + K-P-H
    2  Γ-Σ-M-T’-K-T-Γ-∆-A-R-L-S’-H-S-A
    3  Γ-Σ-M-T’-K-T-Γ-∆-A
    4  Γ-Σ-M
    5  K-T-Γ
sc  1  Γ-∆-X-Y-M-V-R-Λ-Γ-Σ-M
    2  Γ-∆-X-Y-M-V-R-Λ-Γ
    3  Γ-∆-X-Y-M-V-R
    4  Γ-∆-X-Y-M
fcc 1  X-∆-Γ-Λ-L-Q-W-N-K-Σ-Γ
       + L-M-U-S-X-Z-W-D-U
    2  X-∆-Γ-Λ-L-Q-W-N-K-Σ-Γ
    3  X-∆-Γ-Λ-L
    4  Γ-∆-X
    5  Γ-Λ-L
bcc 1  Γ-D-H-G-N-Σ-Γ-Λ-P-F-H + N-D-P
    2  Γ-D-H-G-N-Σ-Γ-Λ-P-F-H
    3  Γ-D-H-G-N-Σ-Γ-Λ-P
    4  Γ-D-H-G-N-Σ-Γ
    5  Γ-D-H
""",
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


def _energy_definition():
    energy = ENERGY(
        emin=(-0.2, "The energy at which the BSF mesh starts", None),
        emax=(-1.0, "The energy at which the BSF mesh ends", None),
        defaults={"GRID": 3, "NE": _ne_default, "ImE": 0.001},
    )
    energy["EMIN"].default_value = _emin_default
    energy["EMAX"].default_value = _emax_default
    energy["EMAXEV"].default_value = _emaxev_default
    return energy


class BSFTaskSection(InputSection):
    def k_path_gui(self, atoms):
        from ase2sprkkr.gui.k_path import k_path_gui

        out = k_path_gui(atoms)
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
