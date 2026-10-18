"""Torque input-parameter definition."""

import numpy as np
from ...common.grammar_types import Integer, SetOf
from ...common.warnings import DataValidityError
from .sections import TASK, CONTROL, TAU, ENERGY, SITES, STRCONST, MODE
from ..input_parameters_definitions import InputParametersDefinition, InputValueDefinition as V
from ..task_definition import TaskDefinition


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



def input_parameters_definition():
    mmode = MODE["MODE"].copy(warning_condition=_torque_mode_warning_condition)
    energy = ENERGY(defaults={"ImE": 0.0, "GRID": 8, "NE": 36})
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = InputParametersDefinition(
        "torque",
        [
            CONTROL("TORQUE"),
            TAU,
            MODE.copy(add=[mmode]),
            STRCONST,
            energy,
            TASK("TORQUE").copy(
                [
                    V("THETAQ", SetOf(float, min_length=1), [90.0], info="Polar angles of the torque axis per site, in degrees."),
                    V("PHIQ", SetOf(float, min_length=1), [90.0], info="Azimuthal angles of the torque axis per site, in degrees."),
                ]
            ),
            SITES,
        ],
        info="TORQUE",
        validators=_validate_torque,
        description="Requires a potential with rotated moments (KMROT=1 or 2), with torque axes perpendicular to them. Native repeats the final angle for remaining sites.",
    )
    schema["CONTROL"]["NONMAG"].warning_condition = _torque_nonmag_warning_condition
    schema._task_definition = TaskDefinition(
        "TORQUE",
        name="TORQUE",
        executable="kkrgen",
        mpi=True,
    )
    return schema
