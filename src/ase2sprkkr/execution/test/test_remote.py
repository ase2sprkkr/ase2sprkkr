import asyncio
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from ..base import ExecutionSpec, FailurePhase, LaunchSpec, RunState, SubmissionNotStarted
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


def test_remote_input_file_and_stdin_staged_only_once(tmp_path):
    transport = _SpyTransport()
    executor, spec, _, remote_root = _text_run(tmp_path, transport)
    runner = executor.create_run(spec)

    async def prepare():
        directory = await runner.executor.map_directory(spec.directory)
        from ..run_model import FileManifest
        return await runner._prepare_request(
            "one-run", directory, FileManifest.from_spec(spec),
        )

    try:
        request = asyncio.run(prepare())
        target = remote_root / "calculation" / "input.txt"
        assert request["input_file"] == str(target)
        assert request["stdin"] == str(target)
        assert request["input_files"] == [str(target)]
        assert [destination for _, destination in transport.put_files] == [target]
    finally:
        spec.input_file.close()
        spec.output_file.close()


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
    "policy, local_value, remote_value, expected, uploads",
    [
        ("changed", "new", "old", "new", 1),
        ("changed", "same", "same", "same", 0),
        ("always", "local", "remote", "local", 1),
        ("missing", "local", "remote", "remote", 0),
        ("never", "local", "local", "local", 0),
    ],
    ids=("changed", "unchanged", "always", "missing", "never"),
)
def test_remote_input_file_policies(
    tmp_path, policy, local_value, remote_value, expected, uploads,
):
    transport = _SpyTransport()
    workdir = tmp_path / "local" / "calculation"
    workdir.mkdir(parents=True)
    artifact = workdir / "artifact.dat"
    artifact.write_text(local_value)
    remote_artifact = tmp_path / "remote" / "calculation" / "artifact.dat"
    remote_artifact.parent.mkdir(parents=True)
    remote_artifact.write_text(remote_value)

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
    assert sum(dst == remote_artifact for _, dst in transport.put_files) == uploads


@pytest.mark.parametrize("direction", ("input", "output"))
def test_remote_never_policy_ignores_missing_declared_file(tmp_path, direction):
    transport = _SpyTransport()
    missing = tmp_path / "local" / "calculation" / "missing.dat"
    kwargs = (
        {"declared_input_files": (missing,), "input_files": "never"}
        if direction == "input" else
        {"declared_output_files": (missing,), "output_files": "never"}
    )
    executor, spec, _, remote_root = _text_run(
        tmp_path, transport, code="print('done')", **kwargs,
    )
    if direction == "input":
        _expose_shared_input(remote_root, "input.txt", b"hello")

    result = _processed_run(executor, spec).run()

    assert result.return_code == 0
    assert not missing.exists()
    if direction == "input":
        assert all(src != missing for src, _ in transport.put_files)
    else:
        assert all(dst != missing for _, dst in transport.get_files)


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


@pytest.mark.parametrize(
    "policy, local_value, generated_value, expected, downloads",
    [
        (FileSyncPolicy.CHANGED, "old", "new", "new", 1),
        (FileSyncPolicy.CHANGED, "same", "same", "same", 0),
        (FileSyncPolicy.ALWAYS, "local", "remote", "remote", 1),
        (FileSyncPolicy.MISSING, "local", "remote", "local", 0),
        (FileSyncPolicy.NEVER, "local", "remote", "local", 0),
    ],
    ids=("changed", "unchanged", "always", "missing", "never"),
)
def test_remote_output_file_policies(
    tmp_path, policy, local_value, generated_value, expected, downloads,
):
    transport = _SpyTransport()
    workdir = tmp_path / "local" / "calculation"
    workdir.mkdir(parents=True)
    artifact = workdir / "artifact.dat"
    artifact.write_text(local_value)

    executor, spec, _, remote_root = _text_run(
        tmp_path,
        transport,
        declared_output_files=(artifact,),
        output_files=policy,
        code=f"import pathlib; pathlib.Path('artifact.dat').write_text({generated_value!r}); print('done')",
    )
    result = _processed_run(executor, spec).run()

    remote_artifact = remote_root / "calculation" / "artifact.dat"
    assert result.return_code == 0
    assert artifact.read_text() == expected
    assert transport.get_files == (
        [(remote_artifact, artifact)] if downloads else []
    )


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


def test_remote_cache_creation_content_and_mtime_changes_in_one_run(tmp_path):
    """Four independent cache outcomes should need only one remote calculation."""
    transport = _SpyTransport()
    workdir = tmp_path / "local" / "calculation"
    workdir.mkdir(parents=True)

    changed = workdir / "changed-cache.dat"
    created = workdir / "created-cache.dat"
    unchanged = workdir / "unchanged-cache.dat"
    touched = workdir / "mtime-cache.dat"
    changed.write_text("old")
    unchanged.write_text("same")
    touched.write_text("same")

    executor, spec, _, remote_root = _text_run(
        tmp_path, transport,
        declared_cache_files=(changed, created, unchanged, touched),
        code=(
            "import os; from pathlib import Path; "
            "changed=Path('changed-cache.dat'); "
            "assert changed.read_text()=='old'; changed.write_text('new'); "
            "Path('created-cache.dat').write_text('generated'); "
            "assert Path('unchanged-cache.dat').read_text()=='same'; "
            "touched=Path('mtime-cache.dat'); "
            "assert touched.read_text()=='same'; "
            "info=touched.stat(); "
            "os.utime(touched, ns=(info.st_atime_ns, info.st_mtime_ns+1000000000)); "
            "print('done')"
        ),
    )
    result = _processed_run(executor, spec).run()
    remote_directory = remote_root / "calculation"

    assert result.return_code == 0
    assert changed.read_text() == "new"
    assert created.read_text() == "generated"
    assert unchanged.read_text() == "same"
    assert touched.read_text() == "same"

    # Existing cache entries are staged once; a newly created one is not.
    for path in (changed, unchanged, touched):
        remote_path = remote_directory / path.name
        assert sum(dst == remote_path for _, dst in transport.put_files) == 1
    assert all(
        dst != remote_directory / created.name for _, dst in transport.put_files
    )

    # Content or mtime changes, and newly created caches, are downloaded.
    for path in (changed, created, touched):
        remote_path = remote_directory / path.name
        assert sum(
            (src, dst) == (remote_path, path)
            for src, dst in transport.get_files
        ) == 1
    assert all(
        src != remote_directory / unchanged.name
        for src, _ in transport.get_files
    )


@pytest.mark.parametrize("direction", ("input", "output"))
def test_remote_cache_never_policy_skips_only_requested_direction(tmp_path, direction):
    transport = _SpyTransport()
    cache = tmp_path / "local" / "calculation" / "cache.dat"
    cache.parent.mkdir(parents=True)
    cache.write_text("local")

    kwargs = {"input_files": "never"} if direction == "input" else {"output_files": "never"}
    executor, spec, _, remote_root = _text_run(
        tmp_path, transport, declared_cache_files=(cache,),
        code="from pathlib import Path; Path('cache.dat').write_text('remote'); print('done')",
        **kwargs,
    )
    if direction == "input":
        _expose_shared_input(remote_root, "input.txt", b"hello")
        _expose_shared_input(remote_root, "cache.dat", "local")
    remote_cache = remote_root / "calculation" / "cache.dat"

    result = _processed_run(executor, spec).run()

    assert result.return_code == 0
    if direction == "input":
        assert all(dst != remote_cache for _, dst in transport.put_files)
        assert cache.read_text() == "remote"
        assert (remote_cache, cache) in transport.get_files
    else:
        assert cache.read_text() == "local"
        assert all(src != remote_cache for src, _ in transport.get_files)


def test_remote_sync_policy_rejects_unknown_value(tmp_path):
    endpoint = _executor()
    with pytest.raises(ValueError):
        endpoint._context(
            endpoint,
            _launch(tmp_path, tmp_path / "remote", input_files="sometimes"),
        )


@pytest.mark.parametrize(
    "root_mode, remote_root, expected",
    [
        ("default", None, (Path.home() / "calculations" / "Fe").as_posix()),
        ("remote-only", "/scratch/tester", "/scratch/tester/calculations/Fe"),
        ("both", "/scratch/tester/projects", "/scratch/tester/projects/project/Fe"),
    ],
)
def test_remote_directory_mapping(tmp_path, root_mode, remote_root, expected):
    endpoint = _executor()
    if root_mode == "both":
        local_root = tmp_path / "local"
        local = local_root / "project" / "Fe"
        launch = LaunchSpec(local_root=local_root, remote_root=remote_root)
    else:
        local = Path.home() / "calculations" / "Fe"
        launch = LaunchSpec() if remote_root is None else LaunchSpec(remote_root=remote_root)
    context = endpoint._context(endpoint, launch)

    mapped = asyncio.run(context.map_directory(local))

    assert mapped == expected


def test_remote_directory_mapping_rejects_path_outside_local_root(tmp_path):
    local_root = tmp_path / "local"
    endpoint = _executor()
    context = endpoint._context(
        endpoint,
        LaunchSpec(local_root=local_root, remote_root="/scratch/tester"),
    )

    with pytest.raises(ValueError, match="outside mapped root"):
        asyncio.run(context.map_directory(tmp_path / "elsewhere"))


def test_remote_none_stream_files_are_live_only(tmp_path):
    workdir = tmp_path / "calculation"
    workdir.mkdir()
    executor = _executor()
    spec = ExecutionSpec(
        argv=(sys.executable, "-c",
              "import sys,time; time.sleep(.1); print('out',flush=True); "
              "print('err',file=sys.stderr,flush=True)"),
        input_file=None,
        output_file=None,
        error_file=None,
        directory=str(workdir),
        launch=_launch(tmp_path, tmp_path / "remote"),
    )
    observed = []
    runner = executor.create_run(spec)
    runner.add_output_listener(lambda text, kind: observed.append((kind, text)))
    outcome = runner.run()
    assert runner.state == RunState.FINISHED
    assert outcome.return_code == 0
    assert outcome.output_file is None
    assert outcome.error_file is None
    assert ("out", "out\n") in observed
    assert ("err", "err\n") in observed
    assert list((tmp_path / "remote").rglob("stdout")) == []
    assert list((tmp_path / "remote").rglob("stderr")) == []


@pytest.mark.parametrize(
    "policy,source,destination,expected",
    [
        (FileSyncPolicy.NEVER, "new", None, False),
        (FileSyncPolicy.NEVER, "new", "old", False),
        (FileSyncPolicy.ALWAYS, "new", "new", True),
        (FileSyncPolicy.ALWAYS, None, "old", False),
        (FileSyncPolicy.MISSING, "new", None, True),
        (FileSyncPolicy.MISSING, "new", "old", False),
        (FileSyncPolicy.CHANGED, "new", None, True),
        (FileSyncPolicy.CHANGED, "new", "old", True),
        (FileSyncPolicy.CHANGED, "new", "new", False),
        (FileSyncPolicy.CHANGED, None, None, False),
        (FileSyncPolicy.CHANGED, None, "old", False),
    ],
)
def test_file_sync_transfer_decision_is_direction_independent(
    policy, source, destination, expected,
):
    assert policy.should_transfer(source, destination) is expected


def test_remote_argv_reuses_staged_path_identity(tmp_path):
    """Relative, absolute and Path arguments share the staging resolver."""
    from ..remote import RemoteRunner

    work = tmp_path / "work"
    source = work / "input.dat"
    runner = object.__new__(RemoteRunner)
    runner.spec = SimpleNamespace(directory=str(work))
    arguments = ("program", "input.dat", source, "--verbose", "./missing.dat")
    staged = {str(source): "/remote/work/input.dat"}

    assert runner._remote_argv(arguments, staged) == [
        "program", "/remote/work/input.dat", "/remote/work/input.dat",
        "--verbose", "./missing.dat",
    ]


@pytest.mark.parametrize("preflight", (False, True))
def test_remote_preparation_failure_has_one_classification(tmp_path, preflight):
    """Direct submission and persisted preflight must report the same failure."""
    directory = tmp_path / "local"
    directory.mkdir()
    spec = ExecutionSpec(
        argv=("ignored",), directory=str(directory),
        input_file=str(directory / "missing.inp"), output_file=None,
        launch=LaunchSpec(
            local_root=str(tmp_path), remote_root=str(tmp_path / "remote"),
            input_sync="always",
        ),
    )
    runner = RemoteEndpointExecutor(LocalTransport()).create_run(spec)

    with pytest.raises(SubmissionNotStarted, match="Local input file"):
        if preflight:
            runner.prepare_submission()
        else:
            runner.submit()

    assert runner.failure_phase == FailurePhase.SUBMISSION
    assert "Local input file" in runner.failure_message
    assert runner.backend_id is None
    assert runner.state == RunState.NEW


@pytest.mark.parametrize("preflight", (False, True))
def test_remote_submit_failure_is_not_misclassified_as_failed_preparation(
    tmp_path, monkeypatch, preflight,
):
    directory = tmp_path / "local"
    directory.mkdir()
    spec = ExecutionSpec(
        argv=("ignored",), directory=str(directory),
        input_file=None, output_file=None,
        launch=LaunchSpec(
            local_root=str(tmp_path), remote_root=str(tmp_path / "remote"),
        ),
    )
    executor = RemoteEndpointExecutor(LocalTransport())
    runner = executor.create_run(spec)

    async def fail_after_preparation(_request):
        raise OSError("backend submission uncertain")

    monkeypatch.setattr(executor, "_submit", fail_after_preparation)
    if preflight:
        runner.prepare_submission()

    with pytest.raises(OSError, match="backend submission uncertain"):
        runner.submit()

    assert runner.failure_phase is None
    assert runner.backend_id is None
