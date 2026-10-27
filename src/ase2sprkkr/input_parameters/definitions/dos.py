"""DOS input-parameter definition."""

from ...common.grammar_types import Integer, SetOf
from ._energy import validate_dos_energy as _validate_dos
from ._sections import CONTROL, TAU, ENERGY, TASK, SITES, MODE, CPA, STRCONST
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import TaskDefinition, TaskFile


def input_parameters_definition():
    energy = ENERGY(
        emax=(1.0, "value of the highest energy", None),
        defaults={"GRID": 3, "NE": 300, "ImE": 0.01},
    )
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = InputParametersDefinition(
        "dos",
        [
            CONTROL("DOS"),
            TAU,
            energy,
            MODE,
            CPA,
            STRCONST,
            TASK("DOS"),
            SITES,
        ],
        info="DOS - The density of states computation",
        validators=_validate_dos,
    )
    schema._task_definition = TaskDefinition(
        "DOS",
        name="DOS",
        executable="kkrgen",
        mpi=True,
        output_reader="DOS",
        output_files=(TaskFile("{dataset}.dos", type="DOS"),),
    )
    return schema
