from ..definitions._file_manifests import tau_file
from ..task_definition import TaskDefinition, TaskFile


def task_definition():
    return TaskDefinition(
        "BSF",
        name="BSFKK",
        executable="kkrgen",
        mpi=True,
        output_reader="BSF",
        cache_files=(TaskFile(tau_file, info="Reusable TAU-directory cache file"),),
    )
