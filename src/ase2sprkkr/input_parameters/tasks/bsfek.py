from ..definitions._file_manifests import bsf_file, tau_file
from ..task_definition import FileSpec, TaskDefinition


def task_definition():
    return TaskDefinition(
        "BSF",
        name="BSFEK",
        executable="kkrgen",
        mpi=True,
        output_reader="BSF",
        output_file_templates=(FileSpec(bsf_file),),
        cache_file_templates=(FileSpec(tau_file),),
    )
