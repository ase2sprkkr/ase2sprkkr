from ..definitions._file_manifests import TauCacheTaskDefinition


def task_definition():
    return TauCacheTaskDefinition(
        "BSF",
        name="BSFKK",
        executable="kkrgen",
        mpi=True,
        result_reader_name="BSF",
    )
