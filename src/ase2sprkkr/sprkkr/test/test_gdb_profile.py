from ase2sprkkr.execution import LaunchSpec
from ase2sprkkr.sprkkr.profile import SPRKKRProfile


def test_gdb_execution_forces_serial_launch(tmp_path):
    profile = SPRKKRProfile(mpi=8, mpi_warning=False)
    input_file = (tmp_path / "run.inp").open("w+")

    try:
        spec = profile.create_execution(
            "kkrscf",
            input_file,
            None,
            str(tmp_path),
            mpi_supported=True,
            gdb=True,
        )
    finally:
        input_file.close()

    assert spec.argv == ("gdb", "kkrscf")
    assert spec.mpi_argv is None
    assert spec.launch == LaunchSpec(
        mpi=False,
        mpi_warning=False,
        input_sync="changed",
        output_sync="changed",
    )
    assert spec.stdin is None
