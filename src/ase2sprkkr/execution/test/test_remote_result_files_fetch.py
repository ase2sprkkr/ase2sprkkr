from ase2sprkkr.execution import LocalTransport, RemoteExecutor, RunReference


def test_remote_executor_fetches_only_requested_historical_files(tmp_path):
    local_root = tmp_path / "local"
    remote_root = tmp_path / "remote"
    local_directory = local_root / "calculation"
    remote_directory = remote_root / "calculation"
    local_directory.mkdir(parents=True)
    remote_directory.mkdir(parents=True)

    remote_input = remote_directory / "run.inp"
    remote_potential = remote_directory / "start.pot"
    remote_output = remote_directory / "run.out"
    remote_first = remote_directory / "result-1.dat"
    remote_second = remote_directory / "result-2.dat"
    remote_later = remote_directory / "result-later.dat"
    remote_input.write_text("input")
    remote_potential.write_text("potential")
    remote_output.write_text("output")
    remote_first.write_text("one")
    remote_second.write_text("two")
    remote_later.write_text("not part of the historical snapshot")

    local_first = local_directory / "result-1.dat"
    local_first.write_text("keep-local")

    executor = RemoteExecutor(LocalTransport())
    reference = RunReference(
        executor="remote:local",
        native_id="opaque-historical-handle",
        directory=str(remote_directory),
        context={
            "launch": {
                "local_root": str(local_root),
                "remote_root": str(remote_root),
                "input_sync": "changed",
                "output_sync": "changed",
            }
        },
    )
    local_input = local_directory / "run.inp"
    local_potential = local_directory / "start.pot"
    local_output = local_directory / "run.out"
    local_second = local_directory / "result-2.dat"
    local_later = local_directory / "result-later.dat"

    assert executor.supports_result_file_retrieval is True
    executor.fetch_result_files(
        reference,
        (local_input, local_potential, local_output, local_second),
    )

    assert local_input.read_text() == "input"
    assert local_potential.read_text() == "potential"
    assert local_output.read_text() == "output"
    assert local_first.read_text() == "keep-local"
    assert local_second.read_text() == "two"
    assert not local_later.exists()
    assert local_output.stat().st_mtime_ns == remote_output.stat().st_mtime_ns
    assert local_second.stat().st_mtime_ns == remote_second.stat().st_mtime_ns
