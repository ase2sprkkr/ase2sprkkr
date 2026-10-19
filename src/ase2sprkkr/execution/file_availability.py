"""Local artifact discovery shared by process and PSI/J runners.

No transport or execution state lives here. The caller decides whether the
producer's output streams are safe to consider published.
"""

from __future__ import annotations

import os

from .transport import expand_local_candidates, has_file_magic


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
        if has_file_magic(path):
            matches = tuple(
                candidate for candidate in expand_local_candidates(path)
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



def validate_output_streams(streams, *, protected=()):
    """Reject a stream that would overwrite a declared input or restart cache."""
    protected_paths = tuple(
        os.path.abspath(os.fspath(p)) for p in protected if p is not None
    )
    for stream in streams:
        if stream is None:
            continue
        target = os.path.abspath(os.fspath(stream))
        for declaration in protected_paths:
            if target in expand_local_candidates(declaration):
                raise ValueError(
                    f"Output stream {target} overlaps protected input/cache {declaration}"
                )


def clear_declared_outputs(paths, *, protected=()):
    """Remove old declared results without destroying inputs or restart caches.

    Called on the execution host before starting a new calculation, never on
    reattach. A wildcard removes only its concrete files.
    """
    protected_paths = (os.path.abspath(os.fspath(p)) for p in protected if p is not None)
    protected = {
        os.path.abspath(candidate)
        for declaration in protected_paths
        for candidate in expand_local_candidates(declaration)
    }
    for declared in dict.fromkeys(os.fspath(p) for p in paths if p is not None):
        for candidate in expand_local_candidates(declared):
            path = os.path.abspath(candidate)
            if path in protected or not os.path.isfile(path):
                continue
            os.unlink(path)
