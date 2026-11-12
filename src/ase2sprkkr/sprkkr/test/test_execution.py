"""SPR-KKR execution declarations and real local KKRunner client runs."""

import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

from kkrunner import ExecutionSpec, RunState
from kkrunner.client import EndpointClient
from kkrunner.execution import MpiSpec
from kkrunner.runs import Run

from ..profile import SPRKKRProfile


def test_profile_builds_serial_command():
    input_file = SimpleNamespace(name="/tmp/Fe.inp")
    output_file = object()
    profile = SPRKKRProfile(mpi=False)
    spec = profile.create_execution(
        "kkrscf", input_file, output_file, "/tmp",
        input_files=(input_file, "/tmp/Fe.pot"),
        output_files=(output_file, "/tmp/Fe_DOS.dos"),
        mpi_supported=True,
    )
    assert spec.argv == ("kkrscf",)
    assert spec.stdin is input_file
    assert spec.input_file is input_file
    assert spec.output_file is output_file
    assert spec.input_files == (input_file, "/tmp/Fe.pot")
    assert spec.output_files == (output_file, "/tmp/Fe_DOS.dos")


def test_profile_builds_portable_mpi_intent():
    input_file = SimpleNamespace(name="/tmp/Fe.inp")
    profile = SPRKKRProfile(
        executable_dir="/opt/sprkkr", executable_suffix="_v97",
        mpi=["mpirun", "-np", "4"],
    )
    spec = profile.create_execution(
        "kkrscf", input_file, object(), "/tmp", mpi_supported=True,
    )
    assert spec.argv == ("/opt/sprkkr/kkrscf_v97",)
    assert spec.mpi_argv == ("/opt/sprkkr/kkrscf_v97MPI", "/tmp/Fe.inp")
    assert spec.mpi == MpiSpec(required=True, launcher=("mpirun", "-np", "4"))
    assert spec.stdin is input_file


def test_non_mpi_task_ignores_profile_mpi():
    input_file = SimpleNamespace(name="/tmp/Fe.inp")
    spec = SPRKKRProfile(mpi=["mpirun", "-np", "4"]).create_execution(
        "kkrgen", input_file, object(), "/tmp", mpi_supported=False,
    )
    assert spec.argv == ("kkrgen",)
    assert spec.stdin is input_file
    assert spec.mpi == MpiSpec()


def _local_run(code):
    directory = tempfile.TemporaryDirectory()
    source = Path(directory.name) / "run.inp"
    source.write_text("input")
    spec = ExecutionSpec(
        argv=(sys.executable, "-S", "-c", code),
        input_file=str(source),
        output_file=str(Path(directory.name) / "run.out"),
        error_file=str(Path(directory.name) / "run.err"),
        directory=directory.name,
    )
    runner = Run.create(EndpointClient(), spec, persistent=False)
    runner._test_directory = directory
    return runner


def test_local_run_streams_output_to_parser_and_listeners():
    runner = _local_run(
        "import sys; print('hello', flush=True); "
        "print('problem', file=sys.stderr, flush=True)"
    )
    seen = []
    runner.add_output_listener(lambda text, kind: seen.append((kind, text)))
    result = runner.run()
    assert runner.state == RunState.FINISHED
    assert Path(result.output_file).read_bytes() == b"hello\n"
    assert Path(result.error_file).read_bytes() == b"problem\n"
    assert result.return_code == 0
    assert ("out", "hello\n") in seen
    assert ("err", "problem\n") in seen


def test_local_run_submit_then_wait():
    runner = _local_run("print('ok')")
    assert runner.state == RunState.NEW
    assert runner.backend_id is None
    runner.submit()
    result = runner.wait()
    assert Path(result.output_file).read_bytes() == b"ok\n"
    assert runner.state == RunState.FINISHED


def test_local_binary_output_uses_replacement_characters():
    runner = _local_run("import os; os.write(1, b'\\xff\\n')")
    seen = []
    runner.add_output_listener(lambda text, kind: seen.append((text, kind)))
    assert runner.run().return_code == 0
    assert seen == [("\ufffd\n", "out")]
