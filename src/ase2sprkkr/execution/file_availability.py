"""Local artifact discovery shared by process and PSI/J runners.

No transport or execution state lives here. The caller decides whether the
producer's output streams are safe to consider published.
"""

from __future__ import annotations

import glob
import os


def discover_local_files(files, *, include_streams=True):
    """Resolve pending local declarations and retain individual wildcard matches."""
    declarations = []
    if include_streams:
        declarations.extend(("output", entry) for entry in (files.stdout, files.stderr))
    declarations.extend(("output", entry) for entry in files.outputs)
    declarations.extend(("cache", item.file) for item in files.caches)

    for kind, entry in declarations:
        if (
            entry is None or entry.local_path is None
            or entry.source_pattern is not None
            or entry.availability_state not in (None, "pending")
        ):
            continue
        path = entry.local_path
        if glob.has_magic(path):
            matches = sorted(
                candidate for candidate in glob.glob(path, recursive=True)
                if os.path.isfile(candidate)
            )
            for concrete in matches:
                files = files.with_file_availability(
                    concrete, "available", kind=kind, source_pattern=path,
                )
            state = "available" if matches else "skipped"
        else:
            state = "available" if os.path.isfile(path) else "skipped"
        files = files.with_file_availability(path, state)
    return files
