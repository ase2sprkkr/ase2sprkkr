"""DMI task using the common JXC/DMI input schema."""

from ..task_definition import TaskDefinition


def task_definition():
    return TaskDefinition("JXC", experimental=True, name="DMI", executable="kkrgen", mpi=True,
        output_reader="JXC",
    )
