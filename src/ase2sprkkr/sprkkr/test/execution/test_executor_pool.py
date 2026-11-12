"""Client pooling keys on transport and backend, not run-specific staging."""

from dataclasses import replace

from kkrunner.client.pool import ClientPool
from kkrunner.files import FileEntry, FileManifest
from kkrunner.runs.conflicts import submission_conflicts
from kkrunner.runs.manager import RunManager
from kkrunner.runs.records import RunStore

from ase2sprkkr.sprkkr.profile import SPRKKRProfile, client_from_profile


def test_unused_local_fields_do_not_change_client_identity(tmp_path):
    plain = SPRKKRProfile()
    unused = SPRKKRProfile(
        host="unused", ssh_options=("-p", "2222"), ssh_zstd=True,
        psij_backend="unused", psij_url="unused",
    )
    assert plain.transport == unused.transport
    assert plain.backend_config == unused.backend_config
    assert client_from_profile(plain) is client_from_profile(unused)
    manager = RunManager(RunStore(tmp_path / "runs"))
    files = FileManifest(stdout=FileEntry(local_path=str(tmp_path / "shared.out")))
    first = manager.new_record(directory=tmp_path, files=files,
                               transport=plain.transport, backend=plain.backend_config)
    second = manager.new_record(directory=tmp_path, files=files,
                                transport=unused.transport, backend=unused.backend_config)
    assert submission_conflicts(first, second)


def test_client_pool_reuses_client_for_distinct_mpi_counts():
    pool = ClientPool()
    first = SPRKKRProfile(mpi=8)
    second = replace(first, mpi=64)
    assert pool.for_config(first.transport, first.backend_config) is pool.for_config(
        second.transport, second.backend_config
    )


def test_ssh_staging_not_part_of_client_identity():
    pool = ClientPool()
    first = SPRKKRProfile(
        executor="ssh", host="cluster", local_root="/home/a",
        remote_root="/scratch/a", input_sync="missing",
    )
    second = replace(
        first, local_root="/home/b", remote_root="/scratch/b", input_sync="always",
    )
    assert first.transport == second.transport
    assert first.backend_config == second.backend_config
    assert first.launch_parameters != second.launch_parameters
    assert pool.for_config(first.transport, first.backend_config) is pool.for_config(
        second.transport, second.backend_config
    )


def test_client_pool_discard_recreates_client():
    pool = ClientPool()
    profile = SPRKKRProfile()
    key = (profile.transport, profile.backend_config)
    first = pool.for_config(*key)
    assert pool.discard(key) is first
    assert pool.for_config(*key) is not first
