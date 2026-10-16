from ..task_definition import TaskDefinition


def task_definition():
    return TaskDefinition(
        "BSF",
        name="BSFEK",
        executable="kkrgen",
        mpi=True,
        output_reader="BSF",
    )
