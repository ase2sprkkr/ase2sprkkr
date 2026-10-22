import asyncio

import pytest

from ..base import SubmissionNotStarted
from ..remote import FileSyncPolicy
from ..transport import LocalTransport
from .test_remote import _processed_run, _text_run


def test_local_transport_put_files_expands_patterns_and_preserves_paths(tmp_path):
    source = tmp_path / "source"
    destination = tmp_path / "destination"
    (source / "a").mkdir(parents=True)
    (source / "b").mkdir()
    (source / "a" / "one.dat").write_text("one")
    (source / "b" / "two.dat").write_text("two")
    (source / "a" / "ignored.txt").write_text("ignored")

    asyncio.run(
        LocalTransport().put_files(
            (source / "a" / "*.dat", source / "b" / "*.dat"),
            source,
            destination,
        )
    )

    assert (destination / "a" / "one.dat").read_text() == "one"
    assert (destination / "b" / "two.dat").read_text() == "two"
    assert not (destination / "a" / "ignored.txt").exists()


def test_local_transport_get_files_expands_source_patterns(tmp_path):
    source = tmp_path / "source"
    destination = tmp_path / "destination"
    (source / "nested").mkdir(parents=True)
    (source / "nested" / "one.out").write_text("one")
    (source / "nested" / "two.out").write_text("two")

    asyncio.run(
        LocalTransport().get_files(
            source / "nested" / "*.out",
            source,
            destination,
        )
    )

    assert (destination / "nested" / "one.out").read_text() == "one"
    assert (destination / "nested" / "two.out").read_text() == "two"


def test_remote_executor_expands_input_wildcards(tmp_path):
    workdir = tmp_path / "local" / "calculation"
    workdir.mkdir(parents=True)
    (workdir / "file-a.pot").write_text("a")
    (workdir / "file-b.pot").write_text("b")
    (workdir / "ignored.txt").write_text("ignored")

    executor, spec, _, remote_root = _text_run(
        tmp_path,
        declared_input_files=(workdir / "file*.pot",),
        code=(
            "import pathlib; "
            "print(':'.join(pathlib.Path(name).read_text() "
            "for name in ('file-a.pot', 'file-b.pot')))"
        ),
    )
    result = _processed_run(executor, spec).run()

    assert result.output == b"a:b\n"
    remote_workdir = remote_root / "calculation"
    assert (remote_workdir / "file-a.pot").read_text() == "a"
    assert (remote_workdir / "file-b.pot").read_text() == "b"
    assert not (remote_workdir / "ignored.txt").exists()


def test_remote_executor_rejects_unmatched_input_wildcard(tmp_path):
    workdir = tmp_path / "local" / "calculation"
    executor, spec, _, _ = _text_run(
        tmp_path,
        declared_input_files=(workdir / "missing*.pot",),
    )

    try:
        _processed_run(executor, spec).run()
    except SubmissionNotStarted as exc:
        assert "input pattern did not match" in str(exc)
    else:
        raise AssertionError("unmatched input wildcard was accepted")


def test_remote_executor_downloads_output_wildcards_and_discards_stale_matches(tmp_path):
    workdir = tmp_path / "local" / "calculation"
    remote_workdir = tmp_path / "remote" / "calculation"
    remote_workdir.mkdir(parents=True)
    (remote_workdir / "artifact-stale.dat").write_text("stale")

    executor, spec, _, _ = _text_run(
        tmp_path,
        declared_output_files=(workdir / "artifact-*.dat",),
        code=(
            "import pathlib; "
            "pathlib.Path('artifact-a.dat').write_text('a'); "
            "pathlib.Path('artifact-b.dat').write_text('b'); print('done')"
        ),
    )
    result = _processed_run(executor, spec).run()

    assert result.return_code == 0
    assert (workdir / "artifact-a.dat").read_text() == "a"
    assert (workdir / "artifact-b.dat").read_text() == "b"
    assert not (workdir / "artifact-stale.dat").exists()
    assert not (remote_workdir / "artifact-stale.dat").exists()


def test_remote_executor_syncs_cache_wildcards_in_both_directions(tmp_path):
    workdir = tmp_path / "local" / "calculation"
    workdir.mkdir(parents=True)
    (workdir / "cache-a.dat").write_text("old")

    executor, spec, _, remote_root = _text_run(
        tmp_path,
        declared_cache_files=(workdir / "cache-*.dat",),
        code=(
            "import pathlib; "
            "a=pathlib.Path('cache-a.dat'); assert a.read_text() == 'old'; "
            "a.write_text('new'); pathlib.Path('cache-b.dat').write_text('created'); "
            "print('done')"
        ),
    )
    result = _processed_run(executor, spec).run()

    assert result.return_code == 0
    assert (workdir / "cache-a.dat").read_text() == "new"
    assert (workdir / "cache-b.dat").read_text() == "created"
    assert (remote_root / "calculation" / "cache-a.dat").read_text() == "new"
    assert (remote_root / "calculation" / "cache-b.dat").read_text() == "created"


def test_exact_and_wildcard_outputs_share_availability_tracking(tmp_path):
    from .test_remote import _SpyTransport

    transport = _SpyTransport()
    workdir = tmp_path / "local" / "calculation"
    exact = workdir / "bands.dat"
    overlap = workdir / "part-a.dat"
    pattern = workdir / "part-*.dat"
    executor, spec, _, remote_root = _text_run(
        tmp_path, transport,
        declared_output_files=(exact, overlap, pattern),
        code=(
            "from pathlib import Path; "
            "Path('bands.dat').write_text('exact'); "
            "Path('part-a.dat').write_text('a'); "
            "Path('part-b.dat').write_text('b'); "
            "Path('part-c.dat').write_text('c')"
        ),
    )
    runner = _processed_run(executor, spec)
    runner.run()

    entries = {entry.local_path: entry for entry in runner.files.outputs}
    assert exact.read_text() == "exact"
    assert overlap.read_text() == "a"
    for name, expected in (("part-b.dat", "b"), ("part-c.dat", "c")):
        path = workdir / name
        assert path.read_text() == expected
        assert entries[str(path)].availability_state == "available"
        assert entries[str(path)].source_pattern == str(pattern)

    assert entries[str(exact)].availability_state == "available"
    assert entries[str(overlap)].availability_state == "available"
    assert entries[str(pattern)].availability_state == "available"

    remote_workdir = remote_root / "calculation"
    for name in ("bands.dat", "part-a.dat", "part-b.dat", "part-c.dat"):
        assert sum(
            src == remote_workdir / name for src, _ in transport.get_files
        ) == 1


def test_exact_and_wildcard_caches_keep_individual_sync_states(tmp_path):
    workdir = tmp_path / "local" / "calculation"
    workdir.mkdir(parents=True)
    exact = workdir / "restart.dat"
    unchanged = workdir / "generated-base.dat"
    pattern = workdir / "generated-*.dat"
    exact.write_text("before")
    unchanged.write_text("unchanged")

    executor, spec, _, _ = _text_run(
        tmp_path,
        declared_cache_files=(exact, pattern),
        output_files=FileSyncPolicy.MISSING,
        code=(
            "from pathlib import Path; "
            "Path('restart.dat').write_text('after'); "
            "Path('generated-new.dat').write_text('new')"
        ),
    )
    runner = _processed_run(executor, spec)
    runner.run()

    entries = {entry.file.local_path: entry.file for entry in runner.files.caches}
    assert exact.read_text() == "before"
    assert unchanged.read_text() == "unchanged"
    assert (workdir / "generated-new.dat").read_text() == "new"
    assert entries[str(exact)].availability_state == "skipped"
    assert entries[str(unchanged)].availability_state == "skipped"
    assert entries[str(pattern)].availability_state == "available"
    new_entry = entries[str(workdir / "generated-new.dat")]
    assert new_entry.availability_state == "available"
    assert new_entry.source_pattern == str(pattern)


def test_remote_executor_tolerates_unmatched_optional_cache_wildcard(tmp_path):
    workdir = tmp_path / "local" / "calculation"
    cache_pattern = workdir / "not-present-*.pot"
    executor, spec, _, remote_root = _text_run(
        tmp_path,
        declared_cache_files=(cache_pattern,),
        code="print('done')",
    )
    runner = _processed_run(executor, spec)
    result = runner.run()

    assert result.return_code == 0
    assert result.output == b"done\n"
    assert not tuple(workdir.glob("not-present-*.pot"))
    assert not tuple((remote_root / "calculation").glob("not-present-*.pot"))


def test_local_transfer_batches_use_their_directional_file_primitives(tmp_path):
    source = tmp_path / "source"
    dest = tmp_path / "dest"
    (source / "nested").mkdir(parents=True)
    (source / "nested" / "result.dat").write_text("data")

    class TracedTransport(LocalTransport):
        def __init__(self):
            self.operations = []

        async def put_file(self, src, dst):
            self.operations.append("put")
            await super().put_file(src, dst)

        async def get_file(self, src, dst):
            self.operations.append("get")
            await super().get_file(src, dst)

    transport = TracedTransport()
    async def check():
        await transport.put_files("nested/*.dat", source, dest)
        assert transport.operations == ["put"]
        transport.operations.clear()
        await transport.get_files("nested/*.dat", source, dest)
        assert transport.operations == ["get", "put"]

    asyncio.run(check())
    assert (dest / "nested" / "result.dat").read_text() == "data"


def test_overlapping_exact_and_wildcard_caches_are_transferred_once(tmp_path):
    from .test_remote import _SpyTransport

    transport = _SpyTransport()
    workdir = tmp_path / "local" / "calculation"
    workdir.mkdir(parents=True)
    exact = workdir / "restart-a.dat"
    exact.write_text("old")
    pattern = workdir / "restart-*.dat"
    executor, spec, _, remote_root = _text_run(
        tmp_path, transport,
        declared_cache_files=(exact, pattern),
        code=(
            "from pathlib import Path; "
            "Path('restart-a.dat').write_text('new'); "
            "Path('restart-b.dat').write_text('created')"
        ),
    )
    runner = _processed_run(executor, spec)
    runner.run()
    remote_a = remote_root / "calculation" / "restart-a.dat"
    remote_b = remote_root / "calculation" / "restart-b.dat"
    assert sum(src == remote_a for src, _ in transport.get_files) == 1
    assert sum(src == remote_b for src, _ in transport.get_files) == 1
    entries = {item.file.local_path: item.file for item in runner.files.caches}
    assert entries[str(exact)].availability_state == "available"
    assert entries[str(pattern)].availability_state == "available"
    assert exact.read_text() == "new"


def test_missing_optional_scf_startup_pattern_is_not_fatal(tmp_path):
    """Successful remote SCF need not emit the optional startup DOS."""
    from ase2sprkkr.input_parameters.task_definition import TaskFile

    workdir = tmp_path / "local" / "calculation"
    executor, spec, _, _ = _text_run(
        tmp_path,
        declared_output_files=(TaskFile(workdir / "*_SCFSTART.dos", optional=True),),
        code="print('done')",
    )
    runner = _processed_run(executor, spec)
    result = runner.run()
    assert result.return_code == 0


@pytest.mark.parametrize("optional_first", [True, False])
def test_overlapping_optional_and_required_wildcards_sync_each_file_once(
    tmp_path, optional_first,
):
    """Distinct wildcard requirements must share one concrete-transfer set."""
    from ase2sprkkr.input_parameters.task_definition import TaskFile
    from .test_remote import _SpyTransport

    transport = _SpyTransport()
    workdir = tmp_path / "local" / "calculation"
    optional = TaskFile(workdir / "part-*.dat", optional=True)
    required = TaskFile(workdir / "part-a*.dat")
    declarations = (optional, required) if optional_first else (required, optional)
    executor, spec, _, remote_root = _text_run(
        tmp_path, transport,
        declared_output_files=declarations,
        output_files=FileSyncPolicy.ALWAYS,
        code=(
            "from pathlib import Path; "
            "Path('part-a.dat').write_text('a'); "
            "Path('part-b.dat').write_text('b')"
        ),
    )
    runner = _processed_run(executor, spec)
    assert runner.run().return_code == 0
    for filename, expected in (("part-a.dat", "a"), ("part-b.dat", "b")):
        assert (workdir / filename).read_text() == expected
        assert sum(
            source == remote_root / "calculation" / filename
            for source, _ in transport.get_files
        ) == 1
    entries = {item.local_path: item for item in runner.files.outputs}
    assert entries[str(optional)].availability_state == "available"
    assert entries[str(required)].availability_state == "available"


def test_missing_optional_exact_remote_output_is_not_fatal(tmp_path):
    """Manifest optionality applies equally to explicit filenames and globs."""
    from ase2sprkkr.input_parameters.task_definition import TaskFile

    workdir = tmp_path / "local" / "calculation"
    optional_file = workdir / "optional-result.dat"
    executor, spec, _, _ = _text_run(
        tmp_path,
        declared_output_files=(TaskFile(optional_file, optional=True),),
        code="print('done')",
    )
    runner = _processed_run(executor, spec)
    result = runner.run()
    assert result.return_code == 0
    entries = {entry.local_path: entry for entry in runner.files.outputs}
    assert entries[str(optional_file)].availability_state == "skipped"


@pytest.mark.parametrize("pattern", [False, True])
def test_remote_missing_optional_input_is_not_staged(tmp_path, pattern):
    """A declared optional input is not sent when neither side has it."""
    from ase2sprkkr.input_parameters.task_definition import TaskFile

    workdir = tmp_path / "local" / "calculation"
    pathname = "optional-*.pot" if pattern else "optional.pot"
    missing = workdir / pathname
    executor, spec, _, remote_root = _text_run(
        tmp_path, declared_input_files=(TaskFile(missing, optional=True),),
        code="print('ok')",
    )
    runner = _processed_run(executor, spec)
    assert runner.run().return_code == 0
    assert not tuple((remote_root / "calculation").glob("optional*.pot"))


def test_mandatory_input_alias_cannot_be_made_optional(tmp_path):
    """A required input wins over an optional declaration of the same file."""
    from ase2sprkkr.input_parameters.task_definition import TaskFile

    workdir = tmp_path / "local" / "calculation"
    missing = workdir / "mandatory.pot"
    executor, spec, _, _ = _text_run(
        tmp_path,
        declared_input_files=(TaskFile(missing, optional=True), missing),
    )
    with pytest.raises(SubmissionNotStarted, match="Conflicting optionality"):
        _processed_run(executor, spec).run()


def test_remote_collects_exact_cache_created_during_run(tmp_path):
    """An absent cache is not an input, but can become a new output."""
    from .test_remote import _SpyTransport

    transport = _SpyTransport()
    workdir = tmp_path / "local" / "calculation"
    cache_file = workdir / "fresh-cache.dat"
    executor, spec, _, remote_root = _text_run(
        tmp_path, transport,
        declared_cache_files=(cache_file,),
        code=(
            "from pathlib import Path; "
            "Path('fresh-cache.dat').write_text('fresh')"
        ),
    )
    runner = _processed_run(executor, spec)
    assert runner.run().return_code == 0
    assert cache_file.read_text() == "fresh"
    remote_cache = remote_root / "calculation" / "fresh-cache.dat"
    assert all(source != cache_file for source, _ in transport.put_files)
    assert sum(source == remote_cache for source, _ in transport.get_files) == 1
    entries = {entry.file.local_path: entry.file for entry in runner.files.caches}
    assert entries[str(cache_file)].availability_state == "available"
