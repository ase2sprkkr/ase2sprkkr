from ase2sprkkr.common.file_versions import FileState, FileVersion
from ase2sprkkr.execution import (
    LocalTransport, RemoteExecutor, RunRecord, RunState, Execution,
    WorkingDirectory, FileManifest,
)
from ase2sprkkr.execution.run_model import FileMapping
from ase2sprkkr.execution.endpoint import endpoint_for_executor


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
    run = RunRecord(
        id="historical", state=RunState.FINISHED,
        created_at="2026-10-13T08:00:00+00:00",
        updated_at="2026-10-13T08:00:00+00:00",
        task="DOS",
        execution=Execution(
            endpoint_for_executor(executor),
            WorkingDirectory(str(local_directory), str(remote_directory)),
            "opaque-historical-handle",
        ),
        files=FileManifest(mapping=FileMapping(str(local_root), str(remote_root))),
    )
    local_input = local_directory / "run.inp"
    local_potential = local_directory / "start.pot"
    local_output = local_directory / "run.out"
    local_second = local_directory / "result-2.dat"
    local_later = local_directory / "result-later.dat"

    assert executor.supports_result_file_retrieval is True
    executor.fetch_result_files(
        run,
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



def test_remote_fetch_uses_transport_expansion_for_exact_and_wildcard(tmp_path):
    import os

    class SpyTransport(LocalTransport):
        def __init__(self):
            self.expanded = []
            self.downloaded = []

        async def expand_files(self, sources, root):
            self.expanded.append(os.fspath(sources))
            return await super().expand_files(sources, root)

        async def get_file(self, source, destination):
            self.downloaded.append(os.fspath(source))
            return await super().get_file(source, destination)

    local_root = tmp_path / "local"
    remote_root = tmp_path / "remote"
    (local_root / "calc").mkdir(parents=True)
    remote_dir = remote_root / "calc"
    remote_dir.mkdir(parents=True)

    (remote_dir / "exact.dat").write_text("exact")
    (remote_dir / "part-a.dat").write_text("a")
    (remote_dir / "part-b.dat").write_text("b")
    (remote_dir / ".hidden.dat").write_text("hidden")
    existing = local_root / "calc" / "part-a.dat"
    existing.write_text("keep")

    transport = SpyTransport()
    executor = RemoteExecutor(transport)
    run = RunRecord(
        id="wildcard-fetch", state=RunState.FINISHED,
        created_at="2026-10-13T08:00:00+00:00",
        updated_at="2026-10-13T08:00:00+00:00",
        task="DOS",
        execution=Execution(
            endpoint_for_executor(executor),
            WorkingDirectory(str(local_root / "calc"), str(remote_dir)),
            "opaque-handle",
        ),
        files=FileManifest(mapping=FileMapping(str(local_root), str(remote_root))),
    )

    executor.fetch_result_files(
        run, (local_root / "calc" / "exact.dat",
              local_root / "calc" / "part-*.dat"),
    )
    assert transport.expanded == [
        str(remote_dir / "exact.dat"), str(remote_dir / "part-*.dat"),
    ]
    assert transport.downloaded == [
        str(remote_dir / "exact.dat"), str(remote_dir / "part-b.dat"),
    ]
    assert (local_root / "calc" / "exact.dat").read_text() == "exact"
    assert existing.read_text() == "keep"
    assert (local_root / "calc" / "part-b.dat").read_text() == "b"
    assert not (local_root / "calc" / ".hidden.dat").exists()


def test_remote_staged_result_keeps_canonical_file_absent(tmp_path):
    local_root = tmp_path / "local"
    remote_root = tmp_path / "remote"
    local_dir = local_root / "calc"
    remote_dir = remote_root / "calc"
    local_dir.mkdir(parents=True)
    remote_dir.mkdir(parents=True)
    local = local_dir / "run.out"
    local.write_bytes(b"original result")
    captured = FileVersion.capture(local, require_digest=True)
    local.unlink()
    (remote_dir / "run.out").write_bytes(b"original result")

    executor = RemoteExecutor(LocalTransport())
    run = RunRecord(
        id="staged-result", state=RunState.FINISHED,
        created_at="2026-10-13T08:00:00+00:00",
        updated_at="2026-10-13T08:00:00+00:00",
        task="DOS",
        execution=Execution(
            endpoint_for_executor(executor),
            WorkingDirectory(str(local_dir), str(remote_dir)),
            "backend-id",
        ),
        files=FileManifest(
            mapping=FileMapping(str(local_root), str(remote_root)),
            output_versions=(captured,),
            captured_at="2026-10-13T08:00:00+00:00",
        ),
    )
    staged = local_dir / ".staged-run.out"
    executor.fetch_result_files(
        run, (local,), destinations={str(local): str(staged)},
    )
    assert not local.exists()
    assert captured.state() is FileState.MISSING
    assert staged.read_bytes() == b"original result"
