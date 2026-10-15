"""File provenance for path-backed values exposed by calculation results."""

from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType
from typing import Any

from ..common.file_versions import FileVersion


def capture_result_file_versions(result: Any) -> Mapping[str, FileVersion]:
    """Capture path-backed output values exposed by *result*.

    Output options are owned by ASE2SPRKKR and expose ``path()`` when their
    value is backed by a file. Provenance enrichment is deliberately best-effort:
    an optional or unstable output is omitted rather than represented by a weak
    metadata-only identity.
    """
    try:
        items = tuple(result.output_values.items())
    except Exception:  # noqa: BLE001 - optional result metadata boundary
        return MappingProxyType({})

    versions: dict[str, FileVersion] = {}
    for key, value in items:
        path_getter = getattr(value, "path", None)
        if not callable(path_getter):
            continue
        try:
            path = path_getter()
        except Exception:  # noqa: BLE001 - optional result metadata boundary
            continue
        if not path:
            continue
        try:
            versions[str(key)] = FileVersion.capture(path, require_digest=True)
        except Exception:  # noqa: BLE001 - provenance enrichment is best effort
            continue
    return MappingProxyType(versions)


__all__ = ["capture_result_file_versions"]
