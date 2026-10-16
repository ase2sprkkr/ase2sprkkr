"""XAS input-parameter definition."""

from ...common.grammar_types import DefKeyword, flag
from ._file_manifests import tau_file
from .sections import CONTROL, TASK, TAU, ENERGY, SITES, MODE
from ..input_parameters_definitions import InputParametersDefinition, InputValueDefinition as V
from ..task_definition import FileSpec, TaskDefinition, TaskFiles


def input_parameters_definition():
    """XAS task input-parameter definition."""
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
                        info="Polar angle θ (FRAMETET) defining the orientation of the electric field vector of the incident light with respect to the material surface normal. Default is FRAMETET = 0, meaning the field lies along the surface normal.",
                    ),
                    V(
                        "FRAMEPHI",
                        float,
                        0.0,
                        is_expert=True,
                        is_optional=True,
                        info="Azimuthal angle φ (FRAMEPHI) defining the in-plane rotation of the electric field vector relative to the surface reference axis. Default is FRAMEPHI = 0, corresponding to alignment with the x-axis of the surface frame.",
                    ),
                ]
            ),
            SITES,
        ],
        info="Calculates X-ray absorption spectra",
    )
    schema._task_definition = TaskDefinition(
        "XAS",
        name="XAS",
        executable="kkrgen",
        mpi=True,
        files=TaskFiles(cache=(FileSpec(tau_file),)),
    )
    return schema
