from ase2sprkkr.input_parameters.task_definition import DeclaredFile
import asyncio
from pathlib import Path

from ..base import SubmissionNotStarted
from ..transport import LocalTransport
from .test_remote import _Reader, _text_run


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
    result = executor.create_run(_Reader(), spec).run()

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
        executor.create_run(_Reader(), spec).run()
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
    result = executor.create_run(_Reader(), spec).run()

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
    result = executor.create_run(_Reader(), spec).run()

    assert result.return_code == 0
    assert (workdir / "cache-a.dat").read_text() == "new"
    assert (workdir / "cache-b.dat").read_text() == "created"
    assert (remote_root / "calculation" / "cache-a.dat").read_text() == "new"
    assert (remote_root / "calculation" / "cache-b.dat").read_text() == "created"


def test_remote_executor_allows_missing_optional_outputs(tmp_path):
    workdir = tmp_path / "local" / "calculation"
    executor, spec, _, _ = _text_run(
        tmp_path,
        declared_output_files=(
            DeclaredFile(workdir / "missing-*.dat", optional=True),
            DeclaredFile(workdir / "not-created.dat", optional=True),
        ),
        code="print('done')",
    )
    result = executor.create_run(_Reader(), spec).run()
    assert result.return_code == 0
    assert not (workdir / "not-created.dat").exists()


def test_remote_executor_allows_missing_optional_input_pattern(tmp_path):
    workdir = tmp_path / "local" / "calculation"
    executor, spec, _, _ = _text_run(
        tmp_path,
        declared_input_files=(DeclaredFile(workdir / "missing-*.pot", optional=True),),
        code="print('done')",
    )
    result = executor.create_run(_Reader(), spec).run()
    assert result.return_code == 0
