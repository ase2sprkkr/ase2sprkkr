"""Public cross-package runtime contracts."""


def test_runtime_contracts_are_available_from_root_api():
    from ase2sprkkr import (
        FileState,
        FileVersion,
        FileVersionCaptureError,
        PotentialFileChanged,
    )
    from ase2sprkkr.common.file_versions import (
        FileState as InternalFileState,
        FileVersion as InternalFileVersion,
        FileVersionCaptureError as InternalFileVersionCaptureError,
    )
    from ase2sprkkr.potentials.potentials import (
        PotentialFileChanged as InternalPotentialFileChanged,
    )

    assert FileState is InternalFileState
    assert FileVersion is InternalFileVersion
    assert FileVersionCaptureError is InternalFileVersionCaptureError
    assert PotentialFileChanged is InternalPotentialFileChanged
