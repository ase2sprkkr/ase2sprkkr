from dataclasses import asdict

import pytest

from ase2sprkkr.execution.endpoint import ExecutorEndpoint
from ase2sprkkr.execution.run_model import (
    CacheEntry, Execution, FileEntry, FileManifest, WorkingDirectory,
)


def test_execution_roundtrip_includes_backend_and_both_directories():
    original = Execution(
        ExecutorEndpoint(executor="ssh", host="cluster"),
        WorkingDirectory("/local/case", "/scratch/case"),
        backend_id="slurm-42",
    )
    assert Execution(**asdict(original)) == original
    assert original.directory.target == "/scratch/case"


def test_local_working_directory_does_not_need_remote_copy():
    run = Execution(ExecutorEndpoint(), WorkingDirectory("/tmp/work"))
    assert run.backend_id is None
    assert run.directory.target == "/tmp/work"


def test_manifest_roundtrip_includes_streams_and_cache_baseline():
    original = FileManifest(
        inputs=(FileEntry("/local/input", "/remote/input"),),
        outputs=(FileEntry("/local/result", "/remote/result"),),
        caches=(CacheEntry(FileEntry("/local/cache", "/remote/cache"), (100, 123)),),
        stdout=FileEntry("/local/out", "/remote/out"),
        stderr=FileEntry("/local/err", "/remote/err"),
    )
    restored = FileManifest(**asdict(original))
    assert restored == original
    assert restored.caches[0].before == (100, 123)


def test_empty_backend_identity_and_pathless_file_are_invalid():
    with pytest.raises(ValueError, match="backend_id"):
        Execution(ExecutorEndpoint(), WorkingDirectory("/tmp"), backend_id="")
    with pytest.raises(ValueError, match="location"):
        FileEntry()
