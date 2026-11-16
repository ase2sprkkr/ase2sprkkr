"""File provenance for path-backed values exposed by calculation results."""

from __future__ import annotations

from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path
from types import MappingProxyType
from typing import Any
import os

from kkrunner.files import FileEntry, FileManifest, FileVersion


def result_file_version(result: Any, path, *, require_digest=False) -> FileVersion | None:
    """Use persisted provenance for tracked results; capture only new imports.

    An absent version in a historical manifest is not permission to adopt the
    bytes currently at that pathname as an output of the old calculation.
    """
    record = getattr(result, "run_record", None)
    if record is None:
        return FileVersion.capture(path, require_digest=require_digest)
    files = record.files
    if files is None or files.captured_at is None:
        return None
    path = os.path.abspath(os.fspath(path))
    return next((version for version in files.all_versions
                 if os.path.abspath(os.fspath(version.path)) == path), None)


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
            version = result_file_version(result, path, require_digest=True)
            if version is not None:
                versions[str(key)] = version
        except Exception:  # noqa: BLE001 - provenance enrichment is best effort
            continue
    return MappingProxyType(versions)



def task_file_declarations(result, *, kinds=("input", "cache", "output")):
    """Return resolved task file declarations, if historical input is usable.

    These declarations are expectations, not evidence that a particular output
    was written by the run. Keep the parser's observed file registry separate.
    """
    try:
        parameters = result.input_parameters
    except Exception:  # noqa: BLE001 - optional historical input parser boundary
        return None, ()
    if parameters is None:
        return None, ()
    from ..input_parameters.input_parameters import InputParameters

    if not isinstance(parameters, InputParameters):
        return None, ()
    # Historical results may retain an original input schema even when the
    # effective TASK value selects another task. That input is useful for
    # inspection, but its values cannot safely resolve the new task's file
    # declarations. Keep observed output files available independently.
    try:
        parameters.task_definition_required()
    except RuntimeError:
        return parameters, ()
    directory = str(result.directory)
    declarations = tuple(
        (kind, value)
        for kind in kinds
        for value in getattr(parameters, f"{kind}_files")(directory)
    )
    return parameters, declarations



def register_declared_result_files(result, *, warn_missing=True):
    """Add existing task-declared outputs to the parsed result file registry.

    The native OUT FILES table is authoritative for observed names and types.
    A task declaration fills gaps but never overrides a file already registered
    by the parser. Absent files are not registered; required ones can warn.
    This must run after the calculation's outputs have been made available
    locally (for persistent/remote runs, after result-file retrieval).
    """
    import glob
    import warnings

    from ..result_files import ResultFile

    _parameters, declarations = task_file_declarations(result, kinds=("output",))
    if not declarations:
        return ()

    directory = Path(result.directory).expanduser().absolute()

    def local_path(value):
        path = Path(os.fspath(value)).expanduser()
        return Path(os.path.abspath(path if path.is_absolute() else directory / path))

    # A path is a stronger identity than a FILES-table key or extension:
    # preserve parser-owned entries and prevent duplicate buttons in the GUI.
    registered = {
        local_path(option.path()): key
        for key, option in result.files.items().items()
    }
    known_formats = ResultFile.definitions
    task = result.task_name.upper()
    registered_now = []
    missing_required = set()

    def result_key(path):
        """Retain the historical accessors for recognized structured outputs."""
        filename = path.name
        if task == "DOS" and filename.lower().endswith(".dos"):
            return "DOS"
        if task in {"BSFEK", "BSFKK"} and filename.lower().endswith(".bsf"):
            if not filename.lower().endswith("_spol.bsf"):
                return "Bloch-SF"
        if task in {"JXC", "DMI"}:
            if filename.endswith(("_XCPLTEN_Jij.dat", "_J_ij.dat")):
                return "jxc"
            if filename.endswith("_XCPLTEN_Dij.dat"):
                return "dij"
            if filename.endswith("_DMIVEC_Dij.dat"):
                return "dmi"
        return filename

    # The ARPES/SPEC result family already owns SPC registration (including
    # multiple polarizations, BAND formats, and header-only SPLEED). Let it
    # decide what is a useful spectrum rather than creating duplicate entries.
    special_spc_tasks = {"ARPES", "AIPES", "BAND", "SPLEED", "CLPES"}

    for _kind, declaration in declarations:
        declared = local_path(declaration)
        pattern = os.fspath(declared)
        paths = (
            sorted(local_path(match) for match in glob.glob(pattern))
            if glob.has_magic(pattern)
            else [declared]
        )
        found = [path for path in paths if path.is_file()]
        if getattr(result, "run_record", None) is not None:
            # Wildcards must not discover additional outputs of a newer run.
            found = [path for path in found if result_file_version(result, path) is not None]
        if not found and not getattr(declaration, "optional", False):
            missing_required.add(pattern)

        hint = getattr(declaration, "type", None)
        file_type = hint.lower() if hint and hint.lower() in known_formats else None
        for path in found:
            if path in registered:
                if file_type and not result.files[registered[path]].file_type:
                    result.files.set_file_type(registered[path], file_type)
                continue
            if task in special_spc_tasks and path.suffix.lower() == ".spc":
                continue

            key = result_key(path)
            if key in result.files:
                # Different directories can contain files with the same name.
                # Keep the filename recognizable without shadowing an entry.
                key = os.path.relpath(path, directory)
                candidate = key
                index = 2
                while candidate in result.files:
                    key = f"{candidate} ({index})"
                    index += 1

            result.files.add_file(key, str(path), file_type)
            registered[path] = key
            registered_now.append(key)

    if warn_missing and missing_required:
        warnings.warn(
            f"{task} result is missing required output file(s): "
            + ", ".join(sorted(missing_required)),
            RuntimeWarning,
            stacklevel=2,
        )
    return tuple(registered_now)


def capture_result_file_manifest(
    result: Any,
    *,
    output_path: str | os.PathLike | None = None,
    potential_path: str | os.PathLike | None = None,
    result_files: Mapping[str, FileVersion] | None = None,
    output_version: FileVersion | None = None,
) -> FileManifest:
    """Snapshot files declared by the task and observed by the parser.

    InputParameters supplies resolved native input, output and cache files.
    The parsed FILES section provides actual filenames and their roles.
    Output wildcards are frozen to matching regular files at import time.
    Missing literal paths remain declarations; unmatched patterns do not
    become fictitious file paths.

    This does not mutate result.files: a declaration does not prove that a
    file was produced by this particular calculation.
    """
    import glob

    from kkrunner.files import CacheEntry

    directory = Path(result.directory).expanduser().absolute()
    source = output_path if output_path is not None else result.output_file
    if source is None:
        raise ValueError("Cannot create a result file manifest without an OUT file")

    def absolute(value):
        path = Path(os.fspath(value)).expanduser()
        return path.absolute() if path.is_absolute() else (directory / path).absolute()

    output = absolute(source)
    # Each path has one semantic role; observed roles take precedence.
    entries: dict[Path, tuple[str, bool, str | None]] = {}

    def add(value, role, *, optional=False, pattern=None, observed=False):
        if value is None or not os.fspath(value):
            return
        path = absolute(value)
        if path == output:
            return
        if path not in entries or observed:
            previous = entries.get(path)
            entries[path] = (
                role,
                previous[1] if previous is not None else optional,
                previous[2] if previous is not None else pattern,
            )

    parameters, declarations = task_file_declarations(result)
    for kind, declaration in declarations:
        path = absolute(declaration)
        optional = getattr(declaration, "optional", False)
        if glob.has_magic(str(path)):
            for match in sorted(glob.glob(str(path))):
                if Path(match).is_file():
                    add(match, kind, optional=optional, pattern=str(path))
        else:
            add(path, kind, optional=optional)

    if parameters is not None:
        # The primary starting potential is a CONTROL value, not necessarily
        # an auxiliary TaskFile declaration or an item in the OUT FILES table.
        potential = parameters.CONTROL.POTFIL()
        if potential:
            add(potential, "input")

    registered = result.files.items()
    # OutputFilesSection.items() returns a dict, unlike dict.items().
    if isinstance(registered, Mapping):
        registered = registered.items()
    for name, option in registered:
        role = str(name).casefold()
        if role == "output":
            continue
        location = option.path() if callable(getattr(option, "path", None)) else option
        if role in {"input", "potential", "start_potential"}:
            kind = "input"
        elif role == "sfn":
            if result.sfn_generated is True:
                kind = "output"
            elif result.sfn_generated is False:
                kind = "input"
            else:
                kind = entries.get(absolute(location), ("input", False, None))[0]
        elif role in {"taudir", "tau_dir"}:
            kind = "cache"
        else:
            # A path which is a declared input/cache stays there even if
            # SPR-KKR also lists it in the OUT FILES section.
            kind = entries.get(absolute(location), ("output", False, None))[0]
        add(location, kind, observed=True)

    add(result.input_file, "input", observed=True)
    if potential_path is not None:
        potential = absolute(potential_path)
        if potential not in entries:
            add(
                potential,
                "output" if str(result.task_name).casefold() == "scf" else "input",
            )

    known = {absolute(v.path): v for v in (result_files or {}).values()}
    for path in known:
        if path != output and path not in entries:
            add(path, "output")

    def selected(kind):
        return tuple(sorted(path for path, item in entries.items() if item[0] == kind))

    def entry(path):
        kind, optional, pattern = entries[path]
        status = None
        message = None
        if kind in {"output", "cache"}:
            if path.is_file():
                status = "available"
            elif kind == "cache" or optional:
                status = "skipped"
            else:
                status = "failed"
                message = f"Required result file is missing: {path}"
        return FileEntry(
            local_path=str(path), optional=optional, source_pattern=pattern,
            availability_state=status, availability_message=message,
        )

    def versions(paths, *, digest):
        return tuple(
            known.get(path) or FileVersion.capture(path, digest=digest)
            for path in paths if path.is_file()
        )

    inputs, outputs, caches = (selected(kind) for kind in ("input", "output", "cache"))
    return FileManifest(
        inputs=tuple(entry(path) for path in inputs),
        outputs=tuple(entry(path) for path in outputs),
        caches=tuple(CacheEntry(entry(path)) for path in caches),
        stdout=FileEntry(
            local_path=str(output),
            availability_state="available" if output.is_file() else "failed",
            availability_message=(
                None if output.is_file() else f"Result OUT file is missing: {output}"
            ),
        ),
        input_versions=versions(inputs, digest=False),
        output_versions=(
            output_version or known.get(output) or FileVersion.capture(output, require_digest=True),
            *versions(outputs, digest=True),
        ),
        captured_at=datetime.now(timezone.utc).isoformat(),
    )


__all__ = ["capture_result_file_versions", "capture_result_file_manifest", "result_file_version"]
