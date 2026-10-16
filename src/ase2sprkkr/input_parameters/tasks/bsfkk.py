from ..definitions._file_manifests import bsf_file, tau_file
from ..task_definition import FileSpec, TaskDefinition, TaskFiles


def task_definition():
    return TaskDefinition(
        "BSF",
        name="BSFKK",
        executable="kkrgen",
        mpi=True,
        output_reader="BSF",
        files=TaskFiles(
            output=(FileSpec(bsf_file),),
            cache=(FileSpec(tau_file),),
        ),
    )
