"""Canonical JXC task using the common JXC/DMI input schema."""

from ..task_definition import TaskDefinition


def task_definition():
    return TaskDefinition(
        "JXC", name="JXC", executable="kkrgen", mpi=True,
        output_reader="JXC",
    )
