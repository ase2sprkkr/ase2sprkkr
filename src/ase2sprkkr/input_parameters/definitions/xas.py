"""XAS task definition."""

from ...common.grammar_types import DefKeyword, flag
from ._file_manifests import TauCacheTaskDefinition
from .sections import CONTROL, TASK, TAU, ENERGY, SITES, MODE
from ..input_parameters_definitions import InputParametersDefinition, InputValueDefinition as V


def task_definition():
    schema = InputParametersDefinition(
        "xas",
        [
            CONTROL("XAS"),
            TAU,
            MODE,
            ENERGY(
                emax=(4.0, "The real part of the highest E-value", None),
                defaults={"ImE": 0.01, "GRID": 6, "NE": 180},
            ),
            TASK("XAS").copy(
                [
                    V("IT", int, 1, is_required=True, info="atom type IT"),
                    V("CL", str, "2P", is_required=True, info="initial core level shell"),
                    V("MECHECK", flag, False, is_optional=True),
                    V(
                        "OUTPUT",
                        DefKeyword(
                            {
                                "MBARN": "output of absorption coefficient µas µ atom = µ Vuc in [Mbarn] ",
                                "SIGMA": "SIGMA : output as absorptive part of optical conductivity σ = µ c/4π in [10E15/s ]",
                            }
                        ),
                        is_expert=True,
                        is_optional=True,
                        info="write extra output",
                    ),
                    V(
                        "FRAMETET",
                        float,
                        0.0,
                        is_expert=True,
                        is_optional=True,
                        info="Polar angle θ defining the incident electric field orientation relative to the surface normal.",
                    ),
                    V(
                        "FRAMEPHI",
                        float,
                        0.0,
                        is_expert=True,
                        is_optional=True,
                        info="Azimuthal angle φ defining the in-plane incident electric field orientation.",
                    ),
                ]
            ),
            SITES,
        ],
        info="Calculates X-ray absorption spectra",
    )
    return TauCacheTaskDefinition(
        schema,
        name="XAS",
        executable="kkrgen",
        mpi=True,
    )
