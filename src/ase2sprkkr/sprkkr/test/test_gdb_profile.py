from ase2sprkkr.execution import LaunchSpec
from ase2sprkkr.sprkkr.profile import SPRKKRProfile


def _create_spec(tmp_path, *, mpi, mpi_supported=True, gdb=True):
    profile = SPRKKRProfile(mpi=mpi, mpi_warning=False)
    input_file = (tmp_path / "run.inp").open("w+")

    try:
        return profile.create_execution(
            "kkrscf",
            input_file,
            None,
            str(tmp_path),
            mpi_supported=mpi_supported,
            gdb=gdb,
        )
    finally:
        input_file.close()


def test_serial_gdb_execution(tmp_path):
    spec = _create_spec(tmp_path, mpi=False)

    assert spec.argv == ("gdb", "kkrscf")
    assert spec.mpi_argv is None
    assert spec.launch == LaunchSpec(
        mpi=False,
        mpi_warning=False,
        input_sync="changed",
        output_sync="changed",
    )
    assert spec.stdin is None


def test_mpi_gdb_execution_uses_xterm_per_rank(tmp_path):
    spec = _create_spec(tmp_path, mpi=4)

    assert spec.argv == ("gdb", "kkrscf")
    assert spec.mpi_argv == (
        "xterm",
        "-e",
        "gdb",
        "--args",
        "kkrscfMPI",
        str(tmp_path / "run.inp"),
    )
    assert spec.launch == LaunchSpec(
        mpi=True,
        process_count=4,
        mpi_warning=False,
        input_sync="changed",
        output_sync="changed",
    )
    assert spec.stdin is None


def test_auto_gdb_keeps_serial_fallback_command(tmp_path):
    spec = _create_spec(tmp_path, mpi="auto")

    assert spec.argv == ("gdb", "kkrscf")
    assert spec.mpi_argv == (
        "xterm",
        "-e",
        "gdb",
        "--args",
        "kkrscfMPI",
        str(tmp_path / "run.inp"),
    )
    assert spec.launch.mpi == "auto"
    assert spec.stdin is None
