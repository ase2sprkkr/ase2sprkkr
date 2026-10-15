from ..task_definition import TaskDefinition


def task_definition():
    return TaskDefinition(
        "ARPES",
        name="BAND",
        executable="kkrspec",
        result_reader_name="SPEC",
    )
