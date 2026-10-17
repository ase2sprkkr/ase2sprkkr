import asyncio
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from ..base import ExecutionSpec, LaunchSpec, RunState, SubmissionNotStarted
from ..postprocessing import PostprocessingRunner
from ..remote import FileSyncPolicy
from ..remote_endpoint import RemoteEndpointExecutor
from ..transport import LocalTransport


class _SpyTransport(LocalTransport):
    def __init__(self):
        self.put_files = []
        self.get_files = []

    async def put_file(self, source, destination):
        self.put_files.append((Path(source), Path(destination)))
        await super().put_file(source, destination)

    async def get_file(self, source, destination):
        self.get_files.append((Path(source), Path(destination)))
        await super().get_file(source, destination)


def _processed_run(executor, spec):
    """Add test-only result processing above the raw remote runner."""

    def process(outcome):
        output = (
            b""
            if outcome.output_file is None
            else Path(outcome.output_file).read_bytes()
        )
        error = (
            b""
            if outcome.error_file is None
            else Path(outcome.error_file).read_bytes()
        )
        return SimpleNamespace(
            return_code=outcome.return_code,
            output=output,
            error=error,
        )

    return PostprocessingRunner(executor.create_run(spec), process)


def _executor(transport=None, **kwargs):
    return RemoteEndpointExecutor(
        transport or LocalTransport(),
        command=(sys.executable, "-m", "ase2sprkkr.tools.ase2sprkkr", "remote"),
        **kwargs,
    )


def _launch(
    local_root,
    remote_root,
    *,
    input_files=FileSyncPolicy.CHANGED,
    output_files=FileSyncPolicy.CHANGED,
):
    return LaunchSpec(
        local_root=local_root,
        remote_root=remote_root,
        input_sync=input_files,
        output_sync=output_files,
    )


def _text_run(
    tmp_path,
    transport=None,
    *,
    declared_input_files=(),
    declared_output_files=(),
    declared_cache_files=(),
    code=None,
    input_files=FileSyncPolicy.CHANGED,
    output_files=FileSyncPolicy.CHANGED,
    **executor_kwargs,
):
    local_root = tmp_path / "local"
    workdir = local_root / "calculation"
    workdir.mkdir(parents=True, exist_ok=True)
    input_file = (workdir / "input.txt").open("w+b")
    input_file.write(b"hello")
    input_file.seek(0)
    output_file = (workdir / "output.txt").open("wb")
    remote_root = tmp_path / "remote"
    executor = _executor(transport=transport, **executor_kwargs)
    spec = ExecutionSpec(
        argv=(sys.executable, "-c", code or "print('ok')"),
        input_file=input_file,
        output_file=output_file,
        directory=str(workdir),
        stdin=input_file,
        input_files=(input_file, *declared_input_files),
        output_files=(output_file, *declared_output_files),
        cache_files=declared_cache_files,
        launch=_launch(
            local_root,
            remote_root,
            input_files=input_files,
            output_files=output_files,
        ),
    )
    return executor, spec, workdir, remote_root


def _expose_shared_input(remote_root, name, data):
    """Simulate an input already visible in the remote view of shared storage."""
    remote = remote_root / "calculation" / name
    remote.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(data, str):
        remote.write_text(data)
    else:
        remote.write_bytes(data)
    return remote


def test_remote_executor_local_transport_streams_output(tmp_path):
    input_file = (tmp_path / "input.txt").open("w+b")
    input_file.write(b"hello")
    input_file.seek(0)
    output_file = (tmp_path / "output.txt").open("wb")

    remote_root = tmp_path / "remote"
    executor = _executor()
    spec = ExecutionSpec(
        argv=(
            sys.executable,
            "-c",
            "import sys; print(sys.stdin.buffer.read().decode().upper())",
        ),
        input_file=input_file,
        output_file=output_file,
        directory=str(tmp_path),
        stdin=input_file,
        launch=_launch(tmp_path, remote_root),
    )

    seen = []
    runner = _processed_run(executor, spec)
    runner.add_output_listener(lambda text, kind: seen.append((text, kind)))
    result = runner.run()

    assert runner.state == RunState.FINISHED
    assert result.return_code == 0
    assert result.output == b"HELLO\n"
    assert result.error == b""
    assert seen == [("HELLO\n", "out")]
    assert runner.backend_id is not None
    assert runner.execution.endpoint is not None
    assert Path(runner.execution.directory.remote) == remote_root


def test_remote_executor_stages_declared_input_files(tmp_path):
    workdir = tmp_path / "calculation"
    workdir.mkdir()
    input_file = (workdir / "input.txt").open("w+b")
    input_file.write(b"hello")
    input_file.seek(0)
    potential_file = workdir / "start.pot"
    potential_file.write_text("potential")
    output_file = (workdir / "output.txt").open("wb")

    remote_root = tmp_path / "remote"
    executor = _executor()
    spec = ExecutionSpec(
        argv=(
            sys.executable,
            "-c",
            (
                "import pathlib, sys; "
                "print(sys.stdin.buffer.read().decode() + ':' + "
                "pathlib.Path('start.pot').read_text())"
            ),
        ),
        input_file=input_file,
        output_file=output_file,
        directory=str(workdir),
        stdin=input_file,
        input_files=(input_file, potential_file),
        launch=_launch(tmp_path, remote_root),
    )

    result = _processed_run(executor, spec).run()

    assert result.return_code == 0
    assert result.output == b"hello:potential\n"
    assert (remote_root / "calculation" / "start.pot").read_text() == "potential"


def test_remote_executor_skips_unchanged_input(tmp_path):
    transport = _SpyTransport()
    workdir = tmp_path / "local" / "calculation"
    workdir.mkdir(parents=True)
    potential = workdir / "start.pot"
    potential.write_text("same")
    remote_potential = tmp_path / "remote" / "calculation" / "start.pot"
    remote_potential.parent.mkdir(parents=True)
    remote_potential.write_text("same")

    executor, spec, _, _ = _text_run(
        tmp_path,
        transport,
        declared_input_files=(potential,),
        code="import pathlib; print(pathlib.Path('start.pot').read_text())",
    )
    result = _processed_run(executor, spec).run()

    assert result.output == b"same\n"
    assert all(dst != remote_potential for _, dst in transport.put_files)


def test_remote_executor_uploads_changed_input(tmp_path):
    transport = _SpyTransport()
    workdir = tmp_path / "local" / "calculation"
    workdir.mkdir(parents=True)
    potential = workdir / "start.pot"
    potential.write_text("new")
    remote_potential = tmp_path / "remote" / "calculation" / "start.pot"
    remote_potential.parent.mkdir(parents=True)
    remote_potential.write_text("old")

    executor, spec, _, _ = _text_run(
        tmp_path,
        transport,
        declared_input_files=(potential,),
        code="import pathlib; print(pathlib.Path('start.pot').read_text())",
    )
    result = _processed_run(executor, spec).run()

    assert result.output == b"new\n"
    assert remote_potential.read_text() == "new"
    assert sum(dst == remote_potential for _, dst in transport.put_files) == 1


def test_remote_executor_uses_remote_only_input(tmp_path):
    transport = _SpyTransport()
    workdir = tmp_path / "local" / "calculation"
    workdir.mkdir(parents=True)
    remote_only = workdir / "remote-only.dat"
    remote_file = tmp_path / "remote" / "calculation" / "remote-only.dat"
    remote_file.parent.mkdir(parents=True)
    remote_file.write_text("remote")

    executor, spec, _, _ = _text_run(
        tmp_path,
        transport,
        declared_input_files=(remote_only,),
        code="import pathlib; print(pathlib.Path('remote-only.dat').read_text())",
    )
    result = _processed_run(executor, spec).run()

    assert result.output == b"remote\n"
    assert not remote_only.exists()
    assert all(dst != remote_file for _, dst in transport.put_files)


@pytest.mark.parametrize(
    "policy, expected, uploaded",
    [
        ("changed", "local", True),
        ("always", "local", True),
        ("missing", "remote", False),
        ("never", "local", False),
    ],
)
def test_remote_input_file_policies(tmp_path, policy, expected, uploaded):
    transport = _SpyTransport()
    workdir = tmp_path / "local" / "calculation"
    workdir.mkdir(parents=True)
    artifact = workdir / "artifact.dat"
    artifact.write_text("local")
    remote_artifact = tmp_path / "remote" / "calculation" / "artifact.dat"
    remote_artifact.parent.mkdir(parents=True)
    remote_artifact.write_text("local" if policy == "never" else "remote")

    executor, spec, _, remote_root = _text_run(
        tmp_path,
        transport,
        declared_input_files=(artifact,),
        input_files=policy,
        code="import pathlib; print(pathlib.Path('artifact.dat').read_text())",
    )
    if policy == "never":
        _expose_shared_input(remote_root, "input.txt", b"hello")
    result = _processed_run(executor, spec).run()

    assert result.output == (expected + "\n").encode()
    assert any(dst == remote_artifact for _, dst in transport.put_files) is uploaded


def test_remote_input_policy_never_only_disables_transfer(tmp_path):
    transport = _SpyTransport()
    missing = tmp_path / "local" / "calculation" / "missing.dat"
    executor, spec, _, remote_root = _text_run(
        tmp_path,
        transport,
        declared_input_files=(missing,),
        input_files="never",
        code="print('done')",
    )
    _expose_shared_input(remote_root, "input.txt", b"hello")

    result = _processed_run(executor, spec).run()

    assert result.return_code == 0
    assert all(src != missing for src, _ in transport.put_files)


def test_remote_input_policy_always_requires_local_file(tmp_path):
    missing = tmp_path / "local" / "calculation" / "missing.dat"
    remote = tmp_path / "remote" / "calculation" / "missing.dat"
    remote.parent.mkdir(parents=True)
    remote.write_text("remote")
    executor, spec, _, _ = _text_run(
        tmp_path, declared_input_files=(missing,), input_files="always"
    )

    with pytest.raises(SubmissionNotStarted, match="Local input file"):
        _processed_run(executor, spec).run()


def test_remote_executor_downloads_changed_output(tmp_path):
    transport = _SpyTransport()
    workdir = tmp_path / "local" / "calculation"
    workdir.mkdir(parents=True)
    artifact = workdir / "artifact.dat"
    artifact.write_text("old")

    executor, spec, _, remote_root = _text_run(
        tmp_path,
        transport,
        declared_output_files=(artifact,),
        code="import pathlib; pathlib.Path('artifact.dat').write_text('new'); print('done')",
    )
    result = _processed_run(executor, spec).run()

    remote_artifact = remote_root / "calculation" / "artifact.dat"
    assert result.return_code == 0
    assert artifact.read_text() == "new"
    assert transport.get_files == [(remote_artifact, artifact)]


def test_remote_executor_skips_unchanged_output(tmp_path):
    transport = _SpyTransport()
    workdir = tmp_path / "local" / "calculation"
    workdir.mkdir(parents=True)
    artifact = workdir / "artifact.dat"
    artifact.write_text("same")

    executor, spec, _, _ = _text_run(
        tmp_path,
        transport,
        declared_output_files=(artifact,),
        code="import pathlib; pathlib.Path('artifact.dat').write_text('same'); print('done')",
    )
    result = _processed_run(executor, spec).run()

    assert result.return_code == 0
    assert artifact.read_text() == "same"
    assert transport.get_files == []


@pytest.mark.parametrize(
    "policy, expected, downloaded",
    [
        (FileSyncPolicy.CHANGED, "remote", True),
        (FileSyncPolicy.ALWAYS, "remote", True),
        (FileSyncPolicy.MISSING, "local", False),
        (FileSyncPolicy.NEVER, "local", False),
    ],
)
def test_remote_output_file_policies(tmp_path, policy, expected, downloaded):
    transport = _SpyTransport()
    workdir = tmp_path / "local" / "calculation"
    workdir.mkdir(parents=True)
    artifact = workdir / "artifact.dat"
    artifact.write_text("local")

    executor, spec, _, remote_root = _text_run(
        tmp_path,
        transport,
        declared_output_files=(artifact,),
        output_files=policy,
        code="import pathlib; pathlib.Path('artifact.dat').write_text('remote'); print('done')",
    )
    _processed_run(executor, spec).run()

    remote_artifact = remote_root / "calculation" / "artifact.dat"
    assert artifact.read_text() == expected
    assert any(src == remote_artifact for src, _ in transport.get_files) is downloaded


def test_remote_output_policy_never_only_disables_transfer(tmp_path):
    missing = tmp_path / "local" / "calculation" / "missing-output.dat"
    executor, spec, _, _ = _text_run(
        tmp_path,
        declared_output_files=(missing,),
        output_files="never",
        code="print('done')",
    )

    result = _processed_run(executor, spec).run()

    assert result.return_code == 0
    assert not missing.exists()


def test_remote_executor_does_not_reuse_stale_output(tmp_path):
    transport = _SpyTransport()
    workdir = tmp_path / "local" / "calculation"
    workdir.mkdir(parents=True)
    artifact = workdir / "artifact.dat"
    remote_artifact = tmp_path / "remote" / "calculation" / "artifact.dat"
    remote_artifact.parent.mkdir(parents=True)
    remote_artifact.write_text("stale")

    executor, spec, _, _ = _text_run(
        tmp_path,
        transport,
        declared_output_files=(artifact,),
        code="print('done')",
    )

    with pytest.raises(FileNotFoundError, match="was not created"):
        _processed_run(executor, spec).run()

    assert not remote_artifact.exists()
    assert not artifact.exists()


def test_remote_cache_is_uploaded_and_changed_cache_is_downloaded(tmp_path):
    transport = _SpyTransport()
    cache = tmp_path / "local" / "calculation" / "cache.dat"
    cache.parent.mkdir(parents=True)
    cache.write_text("old")

    executor, spec, _, remote_root = _text_run(
        tmp_path,
        transport,
        declared_cache_files=(cache,),
        code=(
            "import pathlib; p=pathlib.Path('cache.dat'); "
            "assert p.read_text() == 'old'; p.write_text('new'); print('done')"
        ),
    )
    result = _processed_run(executor, spec).run()

    remote_cache = remote_root / "calculation" / "cache.dat"
    assert result.return_code == 0
    assert cache.read_text() == "new"
    assert any(dst == remote_cache for _, dst in transport.put_files)
    assert (remote_cache, cache) in transport.get_files


def test_remote_cache_created_by_run_is_downloaded(tmp_path):
    transport = _SpyTransport()
    cache = tmp_path / "local" / "calculation" / "cache.dat"

    executor, spec, _, remote_root = _text_run(
        tmp_path,
        transport,
        declared_cache_files=(cache,),
        code="import pathlib; pathlib.Path('cache.dat').write_text('generated'); print('done')",
    )
    result = _processed_run(executor, spec).run()

    remote_cache = remote_root / "calculation" / "cache.dat"
    assert result.return_code == 0
    assert cache.read_text() == "generated"
    assert all(dst != remote_cache for _, dst in transport.put_files)
    assert (remote_cache, cache) in transport.get_files


def test_remote_unchanged_cache_is_not_downloaded(tmp_path):
    transport = _SpyTransport()
    cache = tmp_path / "local" / "calculation" / "cache.dat"
    cache.parent.mkdir(parents=True)
    cache.write_text("same")

    executor, spec, _, remote_root = _text_run(
        tmp_path,
        transport,
        declared_cache_files=(cache,),
        code="import pathlib; assert pathlib.Path('cache.dat').read_text() == 'same'; print('done')",
    )
    _processed_run(executor, spec).run()

    remote_cache = remote_root / "calculation" / "cache.dat"
    assert (remote_cache, cache) not in transport.get_files


def test_remote_cache_change_detection_uses_mtime(tmp_path):
    transport = _SpyTransport()
    cache = tmp_path / "local" / "calculation" / "cache.dat"
    cache.parent.mkdir(parents=True)
    cache.write_text("same")

    executor, spec, _, remote_root = _text_run(
        tmp_path,
        transport,
        declared_cache_files=(cache,),
        code=(
            "import os, pathlib; p=pathlib.Path('cache.dat'); s=p.stat(); "
            "os.utime(p, ns=(s.st_atime_ns, s.st_mtime_ns + 1000000000)); print('done')"
        ),
    )
    _processed_run(executor, spec).run()

    remote_cache = remote_root / "calculation" / "cache.dat"
    assert cache.read_text() == "same"
    assert (remote_cache, cache) in transport.get_files


def test_remote_cache_input_never_skips_upload(tmp_path):
    transport = _SpyTransport()
    cache = tmp_path / "local" / "calculation" / "cache.dat"
    cache.parent.mkdir(parents=True)
    cache.write_text("local")

    executor, spec, _, remote_root = _text_run(
        tmp_path,
        transport,
        declared_cache_files=(cache,),
        input_files="never",
        code="import pathlib; pathlib.Path('cache.dat').write_text('remote'); print('done')",
    )
    _expose_shared_input(remote_root, "input.txt", b"hello")
    remote_cache = _expose_shared_input(remote_root, "cache.dat", "local")
    _processed_run(executor, spec).run()

    assert all(dst != remote_cache for _, dst in transport.put_files)
    assert cache.read_text() == "remote"
    assert (remote_cache, cache) in transport.get_files


def test_remote_cache_output_never_skips_download(tmp_path):
    transport = _SpyTransport()
    cache = tmp_path / "local" / "calculation" / "cache.dat"
    cache.parent.mkdir(parents=True)
    cache.write_text("local")

    executor, spec, _, remote_root = _text_run(
        tmp_path,
        transport,
        declared_cache_files=(cache,),
        output_files="never",
        code="import pathlib; pathlib.Path('cache.dat').write_text('remote'); print('done')",
    )
    _processed_run(executor, spec).run()

    remote_cache = remote_root / "calculation" / "cache.dat"
    assert cache.read_text() == "local"
    assert all(src != remote_cache for src, _ in transport.get_files)


def test_remote_sync_policy_rejects_unknown_value(tmp_path):
    endpoint = _executor()
    with pytest.raises(ValueError):
        endpoint._context(
            endpoint,
            _launch(tmp_path, tmp_path / "remote", input_files="sometimes"),
        )


def test_remote_directory_mapping_defaults_to_home():
    endpoint = _executor()
    context = endpoint._context(endpoint, LaunchSpec())
    local = Path.home() / "calculations" / "Fe"

    mapped = asyncio.run(context.map_directory(local))

    assert mapped == (Path.home() / "calculations" / "Fe").as_posix()


def test_remote_directory_mapping_can_replace_remote_root(tmp_path):
    endpoint = _executor()
    context = endpoint._context(
        endpoint,
        LaunchSpec(remote_root="/scratch/tester"),
    )
    local = Path.home() / "calculations" / "Fe"

    mapped = asyncio.run(context.map_directory(local))

    assert mapped == "/scratch/tester/calculations/Fe"


def test_remote_directory_mapping_can_replace_both_roots(tmp_path):
    local_root = tmp_path / "local"
    local = local_root / "project" / "Fe"
    endpoint = _executor()
    context = endpoint._context(
        endpoint,
        LaunchSpec(
            local_root=local_root,
            remote_root="/scratch/tester/projects",
        ),
    )

    mapped = asyncio.run(context.map_directory(local))

    assert mapped == "/scratch/tester/projects/project/Fe"


def test_remote_directory_mapping_rejects_path_outside_local_root(tmp_path):
    local_root = tmp_path / "local"
    endpoint = _executor()
    context = endpoint._context(
        endpoint,
        LaunchSpec(local_root=local_root, remote_root="/scratch/tester"),
    )

    with pytest.raises(ValueError, match="outside mapped root"):
        asyncio.run(context.map_directory(tmp_path / "elsewhere"))
