from ..task_definition import TaskDefinition, TaskFile


def task_definition():
    return TaskDefinition(
        "BSF",
        name="BSFKK",
        executable="kkrgen",
        mpi=True,
        output_reader="BSF",
        output_files=("{DATASET_ADSI}.bsf", TaskFile("{DATASET_ADSI}_spol.bsf", optional=True)),
    )
