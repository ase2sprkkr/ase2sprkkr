"""Torque input-parameter definition."""

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


def input_parameters_definition():
    mmode = MODE["MODE"].copy(warning_condition=_torque_mode_warning_condition)
    schema = InputParametersDefinition(
        "torque",
        [
            CONTROL("TORQUE"),
            TAU,
            MODE.copy(add=[mmode]),
            STRCONST,
            ENERGY(defaults={"EMIN": -0.2, "ImE": 0.0, "GRID": 8, "NE": 36}),
            TASK("TORQUE").copy(
                [
                    V("THETAQ", [90.0], info="the angles characterizing the orientation of the direction û"),
                    V("PHIQ", [90.0], info="the angles characterizing the orientation of the direction û"),
                ]
            ),
            SITES,
        ],
        info="TORQUE",
    )
    schema["CONTROL"]["NONMAG"].warning_condition = _torque_nonmag_warning_condition
    schema._task_definition = TaskDefinition(
        schema,
        name="TORQUE",
        executable="kkrgen",
        mpi=True,
    )
    return schema
