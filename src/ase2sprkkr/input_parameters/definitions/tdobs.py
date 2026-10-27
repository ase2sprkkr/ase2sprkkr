"""TDOBS native generic DOS/observable path"""

from ...common.grammar_types import Integer, SetOf
from ._energy import validate_dos_energy as _validate_tdobs
from ._file_manifests import TauCacheInputParametersDefinition
from ._sections import CONTROL, TAU, ENERGY, TASK, SITES, MODE, CPA, STRCONST
from ..input_parameters_definitions import InputValueDefinition as V
from ..task_definition import TaskDefinition


def input_parameters_definition():
    energy = ENERGY(emax=(1.0, "value of the highest energy", None), defaults={"GRID": 3, "NE": 300, "ImE": 0.01})
    energy.copy_member("NE", type=SetOf(Integer(min=1), min_length=1, max_length=2))
    schema = TauCacheInputParametersDefinition(
        "tdobs",
        [
            CONTROL("TDOBS"),
            TAU.copy(add=[V("TAUDIR", str, is_required=False, info="Stored tau file; native default DATASET.tau.dir."),
                          V("TAUPATH", str, is_required=False, info="Directory containing TAUDIR; native default current directory.")]),
            energy,
            MODE,
            CPA,
            STRCONST,
            TASK("TDOBS"),
            SITES,
        ],
        validators=_validate_tdobs,
        info="TDOBS: native generic energy-loop DOS and observable calculation with attempted reuse of stored tau. "
             "SPR-KKR 9.7 has no separate TDOBS postprocessing branch. A missing tau file triggers its normal "
             "recomputation fallback. Serial only: native GEN MPI guard does not permit TDOBS."
    )
    schema._task_definition = TaskDefinition("TDOBS", experimental=True, executable="kkrgen", mpi=False, output_files=("{dataset}.dos",))
    return schema


""" DOS - density of states input parameters definition"""
