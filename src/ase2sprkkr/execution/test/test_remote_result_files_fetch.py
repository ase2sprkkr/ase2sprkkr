import json

from ase2sprkkr.execution import LocalTransport, RemoteExecutor, RunReference


def test_remote_executor_fetches_missing_historical_result_files(tmp_path):
    local_root = tmp_path / "local"
    remote_root = tmp_path / "remote"
    local_directory = local_root / "calculation"
    remote_directory = remote_root / "calculation"
    local_directory.mkdir(parents=True)
    remote_directory.mkdir(parents=True)

    remote_input = remote_directory / "run.inp"
    remote_potential = remote_directory / "start.pot"
    remote_output = remote_directory / "run.out"
    remote_input.write_text("input")
    remote_potential.write_text("potential")
    remote_output.write_text("output")

    run_id = "historical-run"
    metadata = remote_directory / ".ase2sprkkr" / "runs" / run_id
    metadata.mkdir(parents=True)
    (metadata / "request.json").write_text(
        json.dumps(
            {
                "launch": {
                    "local_root": str(local_root),
                    "remote_root": str(remote_root),
                    "input_sync": "changed",
                    "output_sync": "changed",
                },
                "input_file": str(remote_input),
                "client_sync": {
                    "input_files": [str(remote_input), str(remote_potential)],
                    "output_file": str(remote_output),
                    "output_files": [],
                    "output_patterns": [],
                },
            }
        )
    )

    executor = RemoteExecutor(LocalTransport())
    reference = RunReference(
        executor="remote:local",
        native_id=run_id,
        directory=str(remote_directory),
    )

    executor.fetch_result_files(reference, result_files=None)

    assert (local_directory / "run.inp").read_text() == "input"
    assert (local_directory / "start.pot").read_text() == "potential"
    assert (local_directory / "run.out").read_text() == "output"
