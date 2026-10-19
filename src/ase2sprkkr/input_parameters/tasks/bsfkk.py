from ..task_definition import TaskDefinition, TaskFile


def task_definition():
    return TaskDefinition(
        "BSF",
        name="BSFKK",
        executable="kkrgen",
        mpi=True,
        output_reader="BSF",
        output_files=("{dataset}_BLOCHSF.bsf", TaskFile("{dataset}_BLOCHSF_spol.bsf", optional=True)),
    )
