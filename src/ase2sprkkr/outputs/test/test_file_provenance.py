from types import SimpleNamespace

from ase2sprkkr.common.file_versions import FileState
from ase2sprkkr.outputs import capture_result_file_versions


class _OutputValue:
    def __init__(self, path):
        self._path = path

    def path(self):
        return self._path


def test_capture_result_file_versions_uses_content_identity(tmp_path):
    path = tmp_path / "result.dat"
    path.write_bytes(b"result")
    result = SimpleNamespace(output_values={"DATA": _OutputValue(path)})

    versions = capture_result_file_versions(result)

    version = versions["DATA"]
    assert version.path == path.absolute()
    assert version.digest is not None
    assert version.available
    assert version.state() is FileState.CURRENT


def test_capture_result_file_versions_ignores_non_path_outputs(tmp_path):
    result = SimpleNamespace(output_values={"VALUE": object()})

    assert dict(capture_result_file_versions(result)) == {}
