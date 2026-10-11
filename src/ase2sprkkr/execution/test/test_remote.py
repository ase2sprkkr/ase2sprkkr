import asyncio
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from ..base import ExecutionSpec, RunState
from ..remote import FileSyncPolicy, RemoteExecutor
from ..transport import LocalTransport


class _Reader:
    print_output = False

    def create_result(self, output_file, input_file):
        return SimpleNamespace(return_code=None, output=None, error=None)

    async def parse_streams(self, stdout, stderr, result):
        return await asyncio.gather(stdout.read(), stderr.read())

    def complete_result(self, result, output, error, return_code):
        result.output = output
        result.error = error
        result.return_code = return_code
        return result


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


def _executor(local_root, remote_root, transport=None, **kwargs):
    return RemoteExecutor(
        transport or LocalTransport(),
        local_root=local_root,
        remote_root=remote_root,
        command=(sys.executable, "-m", "ase2sprkkr.tools.ase2sprkkr", "remote"),
        **kwargs,
    )


def _text_run(
    tmp_path,
    transport=None,
    *,
    input_files=(),
    output_files=(),
    code=None,
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
    executor = _executor(
        local_root, remote_root, transport=transport, **executor_kwargs
    )
    spec = ExecutionSpec(
        argv=(sys.executable, "-c", code or "print('ok')"),
        input_file=input_file,
        output_file=output_file,
        directory=str(workdir),
        stdin=input_file,
        input_files=(input_file, *input_files),
        output_files=(output_file, *output_files),
    )
    return executor, spec, workdir, remote_root


def test_remote_executor_local_transport_streams_output(tmp_path):
    input_file = (tmp_path / "input.txt").open("w+b")
    input_file.write(b"hello")
    input_file.seek(0)
    output_file = (tmp_path / "output.txt").open("wb")

    remote_root = tmp_path / "remote"
    executor = _executor(tmp_path, remote_root)
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
    )

    seen = []
    runner = executor.create_run(_Reader(), spec)
    runner.add_output_listener(lambda text, kind: seen.append((text, kind)))
    result = runner.run()

    assert runner.state == RunState.FINISHED
    assert result.return_code == 0
    assert result.output == b"HELLO\n"
    assert result.error == b""
    assert seen == [("HELLO\n", "out")]
    assert runner.reference is not None
    assert runner.reference.executor == "remote:local"
    assert Path(runner.reference.directory) == remote_root


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
    executor = _executor(tmp_path, remote_root)
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
    )

    result = executor.create_run(_Reader(), spec).run()

    assert result.return_code == 0
    assert result.output == b"hello:potential\n"
    assert (remote_root / "calculation" / "start.pot").read_text() == "potential"


def test_remote_executor_skips_unchanged_input(tmp_path):
    transport = _SpyTransport()
    local_root = tmp_path / "local"
    workdir = local_root / "calculation"
    workdir.mkdir(parents=True)
    potential = workdir / "start.pot"
    potential.write_text("same")

    remote_root = tmp_path / "remote"
    remote_workdir = remote_root / "calculation"
    remote_workdir.mkdir(parents=True)
    remote_potential = remote_workdir / "start.pot"
    remote_potential.write_text("same")

    executor, spec, _, _ = _text_run(
        tmp_path,
        transport,
        input_files=(potential,),
        code="import pathlib; print(pathlib.Path('start.pot').read_text())",
    )
    result = executor.create_run(_Reader(), spec).run()

    assert result.output == b"same\n"
    assert all(destination != remote_potential for _, destination in transport.put_files)


def test_remote_executor_uploads_changed_input(tmp_path):
    transport = _SpyTransport()
    local_root = tmp_path / "local"
    workdir = local_root / "calculation"
    workdir.mkdir(parents=True)
    potential = workdir / "start.pot"
    potential.write_text("new")

    remote_root = tmp_path / "remote"
    remote_workdir = remote_root / "calculation"
    remote_workdir.mkdir(parents=True)
    remote_potential = remote_workdir / "start.pot"
    remote_potential.write_text("old")

    executor, spec, _, _ = _text_run(
        tmp_path,
        transport,
        input_files=(potential,),
        code="import pathlib; print(pathlib.Path('start.pot').read_text())",
    )
    result = executor.create_run(_Reader(), spec).run()

    assert result.output == b"new\n"
    assert remote_potential.read_text() == "new"
    assert sum(destination == remote_potential for _, destination in transport.put_files) == 1


def test_remote_executor_uses_remote_only_input(tmp_path):
    transport = _SpyTransport()
    local_root = tmp_path / "local"
    workdir = local_root / "calculation"
    workdir.mkdir(parents=True)
    remote_only = workdir / "remote-only.dat"

    remote_root = tmp_path / "remote"
    remote_workdir = remote_root / "calculation"
    remote_workdir.mkdir(parents=True)
    remote_file = remote_workdir / "remote-only.dat"
    remote_file.write_text("remote")

    executor, spec, _, _ = _text_run(
        tmp_path,
        transport,
        input_files=(remote_only,),
        code="import pathlib; print(pathlib.Path('remote-only.dat').read_text())",
    )
    result = executor.create_run(_Reader(), spec).run()

    assert result.output == b"remote\n"
    assert not remote_only.exists()
    assert all(destination != remote_file for _, destination in transport.put_files)


@pytest.mark.parametrize(
    "policy, expected, uploaded",
    [
        ("changed", "local", True),
        ("always", "local", True),
        ("missing", "remote", False),
        ("never", "remote", False),
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
    remote_artifact.write_text("remote")

    executor, spec, _, _ = _text_run(
        tmp_path,
        transport,
        input_files=(artifact,),
        input_files=policy,
        code="import pathlib; print(pathlib.Path('artifact.dat').read_text())",
    )
    result = executor.create_run(_Reader(), spec).run()

    assert result.output == (expected + "\n").encode()
    assert any(dst == remote_artifact for _, dst in transport.put_files) is uploaded


def test_remote_input_policy_always_requires_local_file(tmp_path):
    missing = tmp_path / "local" / "calculation" / "missing.dat"
    remote = tmp_path / "remote" / "calculation" / "missing.dat"
    remote.parent.mkdir(parents=True)
    remote.write_text("remote")
    executor, spec, _, _ = _text_run(
        tmp_path, input_files=(missing,), input_files="always"
    )

    with pytest.raises(FileNotFoundError, match="Local input file"):
        executor.create_run(_Reader(), spec).run()


def test_remote_executor_downloads_changed_output(tmp_path):
    transport = _SpyTransport()
    local_root = tmp_path / "local"
    workdir = local_root / "calculation"
    workdir.mkdir(parents=True)
    artifact = workdir / "artifact.dat"
    artifact.write_text("old")

    executor, spec, _, remote_root = _text_run(
        tmp_path,
        transport,
        output_files=(artifact,),
        code="import pathlib; pathlib.Path('artifact.dat').write_text('new'); print('done')",
    )
    result = executor.create_run(_Reader(), spec).run()

    remote_artifact = remote_root / "calculation" / "artifact.dat"
    assert result.return_code == 0
    assert artifact.read_text() == "new"
    assert transport.get_files == [(remote_artifact, artifact)]


def test_remote_executor_skips_unchanged_output(tmp_path):
    transport = _SpyTransport()
    local_root = tmp_path / "local"
    workdir = local_root / "calculation"
    workdir.mkdir(parents=True)
    artifact = workdir / "artifact.dat"
    artifact.write_text("same")

    executor, spec, _, _ = _text_run(
        tmp_path,
        transport,
        output_files=(artifact,),
        code="import pathlib; pathlib.Path('artifact.dat').write_text('same'); print('done')",
    )
    result = executor.create_run(_Reader(), spec).run()

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
        output_files=(artifact,),
        output_files=policy,
        code="import pathlib; pathlib.Path('artifact.dat').write_text('remote'); print('done')",
    )
    executor.create_run(_Reader(), spec).run()

    remote_artifact = remote_root / "calculation" / "artifact.dat"
    assert artifact.read_text() == expected
    assert any(src == remote_artifact for src, _ in transport.get_files) is downloaded


def test_remote_executor_does_not_reuse_stale_output(tmp_path):
    transport = _SpyTransport()
    local_root = tmp_path / "local"
    workdir = local_root / "calculation"
    workdir.mkdir(parents=True)
    artifact = workdir / "artifact.dat"

    remote_root = tmp_path / "remote"
    remote_workdir = remote_root / "calculation"
    remote_workdir.mkdir(parents=True)
    remote_artifact = remote_workdir / "artifact.dat"
    remote_artifact.write_text("stale")

    executor, spec, _, _ = _text_run(
        tmp_path,
        transport,
        output_files=(artifact,),
        code="print('done')",
    )

    with pytest.raises(FileNotFoundError, match="was not created"):
        executor.create_run(_Reader(), spec).run()

    assert not remote_artifact.exists()
    assert not artifact.exists()


def test_remote_sync_policy_rejects_unknown_value(tmp_path):
    with pytest.raises(ValueError):
        _executor(tmp_path, tmp_path / "remote", input_files="sometimes")


def test_remote_directory_mapping_defaults_to_home():
    executor = RemoteExecutor(LocalTransport())
    local = Path.home() / "calculations" / "Fe"

    mapped = asyncio.run(executor.map_directory(local))

    assert mapped == (Path.home() / "calculations" / "Fe").as_posix()


def test_remote_directory_mapping_can_replace_remote_root(tmp_path):
    executor = RemoteExecutor(
        LocalTransport(), remote_root="/scratch/tester"
    )
    local = Path.home() / "calculations" / "Fe"

    mapped = asyncio.run(executor.map_directory(local))

    assert mapped == "/scratch/tester/calculations/Fe"


def test_remote_directory_mapping_can_replace_both_roots(tmp_path):
    local_root = tmp_path / "local"
    local = local_root / "project" / "Fe"
    executor = RemoteExecutor(
        LocalTransport(),
        local_root=local_root,
        remote_root="/scratch/tester/projects",
    )

    mapped = asyncio.run(executor.map_directory(local))

    assert mapped == "/scratch/tester/projects/project/Fe"


def test_remote_directory_mapping_rejects_path_outside_local_root(tmp_path):
    local_root = tmp_path / "local"
    executor = RemoteExecutor(
        LocalTransport(),
        local_root=local_root,
        remote_root="/scratch/tester",
    )

    with pytest.raises(ValueError, match="outside mapped root"):
        asyncio.run(executor.map_directory(tmp_path / "elsewhere"))
