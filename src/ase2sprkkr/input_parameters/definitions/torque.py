"""Torque task input parameters definition"""

import numpy as np
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError

from ._sections import TASK, CONTROL, TAU, ENERGY, SITES, STRCONST, MODE
from ..input_parameters_definitions import InputParametersDefinition as InputParameters, InputValueDefinition as V


def _torque_mode_warning_condition(value):
    return (
        "Torque task does not support SREL (scalar relativity without spin) or "
        "NREL (no relativity at all) MODE. Please change SCF.MODE, or "
        "the computation will fail."
        if value in ("SREL", "NREL")
        else None
    )


def _torque_nonmag_warning_condition(value):
    return (
        "Torque task does not support non-magnetic computation. Please disable "
        "CONTROL.NONMAG, or the computation will fail."
        if value
        else None
    )


def _validate_torque(parameters, _values, why):
    if why == "set":
        return None
    for name in ("THETAQ", "PHIQ"):
        if not np.isfinite(parameters.TASK[name]()).all():
            return DataValidityError("TORQUE axis angles must be finite.")
    if parameters.ENERGY.SPLITSS() or len(parameters.ENERGY.GRID()) != 1 or len(parameters.ENERGY.NE()) != 1:
        return DataValidityError("TORQUE uses only the first energy path; split contours are not supported.")
    return None


def input_parameters():
    """Torque -Torque task input parameters definition"""
    mmode = MODE["MODE"].copy(warning_condition=_torque_mode_warning_condition)
    energy = ENERGY(defaults={"ImE": 0.0, "GRID": 8, "NE": 36})
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    input_parameters = InputParameters(
        "torque",
        [
            CONTROL("TORQUE"),
            TAU,
            MODE.copy(add = [ mmode ]),
            STRCONST,
            energy,
            TASK("TORQUE").copy(
                [
                    V("THETAQ", SetOf(float, min_length=1), [90.0], info="Polar angles of the torque axis u in degrees, per site; u must be perpendicular to magnetization."),
                    V("PHIQ", SetOf(float, min_length=1), [90.0], info="Azimuthal angles of the torque axis u in degrees, per site; these do not rotate magnetization."),
                ]
            ),
            SITES,
        ],
        executable="kkrgen",
        mpi=True,
        info="TORQUE",
        validators=_validate_torque,
        description="Requires a potential with rotated moments (KMROT=1 or 2) and torque axes perpendicular to them. "
                    "The native reader repeats the last supplied angle for remaining sites; arrays are indexed by site, not type.",
    )
    input_parameters["CONTROL"]["NONMAG"].warning_condition=_torque_nonmag_warning_condition
    return input_parameters
