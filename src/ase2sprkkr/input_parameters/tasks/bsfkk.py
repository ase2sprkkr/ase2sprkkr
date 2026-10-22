from ..task_definition import TaskDefinition, DeclaredFile


def task_definition():
    return TaskDefinition(
        "BSF",
        name="BSFKK",
        executable="kkrgen",
        mpi=True,
        output_reader="BSF",
        output_file_templates=("{dataset}_{adsi}.bsf", DeclaredFile("{dataset}_{adsi}_spol.bsf", optional=True)),
    )
