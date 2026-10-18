from ..task_definition import TaskDefinition


def task_definition():
    return TaskDefinition(
        "JXC", name="DMI", executable="kkrgen", mpi=True,
        output_reader="JXC",
    )
