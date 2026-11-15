"""Safe persistence of the virtual ASE SPR-KKR execution profile.

ASE configuration is INI, not ASE2SPRKKR YAML. Only the [sprkkr]
section in the highest-priority permitted file defining [sprkkr] is edited;
Linux /etc and Windows shared configuration are protected by default.
"""

from __future__ import annotations

import configparser
import os
import ntpath
from pathlib import Path, PureWindowsPath
import re
import stat
import sys
import tempfile


_SECTION = re.compile(r"^\s*\[([^\]]+)\]\s*(?:[#;].*)?$")
_OPTION = re.compile(r"^(\s*)([^\s=:#;][^=:]*?)(\s*[=:]\s*)(.*)$")


def ase_config_candidates() -> tuple[Path, ...]:
    """Return ASE's configuration search paths in increasing precedence.

    The last listed file overrides earlier ones when ASE reads config.
    An explicitly empty ASE_CONFIG_PATH disables ASE's configuration,
    so it must not silently create a config file at an unrelated path.
    """
    paths = os.environ.get("ASE_CONFIG_PATH")
    if paths is not None:
        separator = ";" if sys.platform == "win32" else ":"
        names = [name for name in paths.split(separator) if name]
        if not names:
            raise ValueError(
                "ASE_CONFIG_PATH is empty; no ASE configuration file is selected"
            )
        return tuple(Path(name).expanduser() for name in names)

    from ase.config import ASE_CONFIG_FILE
    return (Path(ASE_CONFIG_FILE).expanduser(),)


def _has_sprkkr_section(path: Path) -> bool:
    """Check whether one ASE configuration file defines [sprkkr].

    Each file is inspected independently: ASE merges these sections across
    files, so a high-priority file with an unrelated calculator must not
    move edits away from the file already defining SPR-KKR.
    """
    parser = configparser.ConfigParser(interpolation=None)
    with path.open(encoding="utf-8") as stream:
        parser.read_file(stream)
    return parser.has_section("sprkkr")


_SYSTEM_CONFIG_DIRECTORY = Path("/etc")


def _is_system_config_path(path: Path) -> bool:
    """Classify protected global configuration paths, including symlinks."""
    if sys.platform.startswith("linux"):
        root = _SYSTEM_CONFIG_DIRECTORY.resolve()
        return (
            path.absolute().is_relative_to(root)
            or path.resolve().is_relative_to(root)
        )
    if sys.platform == "win32":
        # ProgramData is shared configuration; Users/Default is a template
        # for newly created accounts, *not* this user's private directory.
        candidate = ntpath.normcase(str(PureWindowsPath(str(path))))
        roots = (
            os.environ.get("ProgramData", "C:\\ProgramData"),
            os.environ.get("ALLUSERSPROFILE", "C:\\ProgramData"),
            os.environ.get("SystemRoot", "C:\\Windows"),
            os.environ.get("ProgramFiles", "C:\\Program Files"),
            os.environ.get("ProgramFiles(x86)", "C:\\Program Files (x86)"),
            str(PureWindowsPath(os.environ.get("SystemDrive", "C:") + chr(92)) / "Users" / "Default"),
        )
        for root in roots:
            try:
                if ntpath.commonpath((candidate, ntpath.normcase(root))) == ntpath.normcase(root):
                    return True
            except ValueError:
                continue
        # Check real location as well, to avoid links into system directories.
        resolved = path.resolve()
        if str(resolved) != str(path):
            return _is_system_config_path_resolved(resolved, roots)
    return False


def _is_system_config_path_resolved(path: Path, roots: tuple[str, ...]) -> bool:
    """Check a resolved Windows path without following links again."""
    candidate = ntpath.normcase(str(PureWindowsPath(str(path))))
    for root in roots:
        try:
            normalized = ntpath.normcase(root)
            if ntpath.commonpath((candidate, normalized)) == normalized:
                return True
        except ValueError:
            continue
    return False


def ase_config_target(*, system: bool = False) -> Path:
    """Select the highest-priority *permitted* ASE file with [sprkkr].

    By default, Linux /etc and Windows shared system configuration are
    excluded. Do not write a lower-priority user file if a protected
    higher-priority [sprkkr] section would override the edit.
    """
    paths = ase_config_candidates()
    permitted = [
        (index, path) for index, path in enumerate(paths)
        if system or not _is_system_config_path(path)
    ]
    if not permitted:
        raise PermissionError(
            "No user-level ASE configuration path is available. "
            "Configure a user path in ASE_CONFIG_PATH, or explicitly "
            "request --system (system=True in Python)."
        )

    target_index, target = next(
        ((index, path) for index, path in reversed(permitted)
         if path.is_file() and _has_sprkkr_section(path)),
        permitted[-1],
    )

    if not system:
        for path in paths[target_index + 1:]:
            if (
                _is_system_config_path(path)
                and path.is_file()
                and _has_sprkkr_section(path)
            ):
                raise PermissionError(
                    "A higher-priority system ASE file defines [sprkkr], "
                    "so this user-level edit would be overridden. "
                    "Change ASE_CONFIG_PATH to put a user file last, "
                    "or explicitly request --system."
                )
    return target


def _rewrite_section(original: str, updates: dict[str, str | None]) -> str:
    """Edit INI options without rewriting unrelated sections or comments."""
    lines = original.splitlines(keepends=True)
    first = None
    last = len(lines)
    for i, line in enumerate(lines):
        found = _SECTION.match(line.rstrip("\r\n"))
        if not found:
            continue
        if first is not None:
            last = i
            break
        if found.group(1).strip().lower() == "sprkkr":
            first = i

    if first is None:
        if not any(value is not None for value in updates.values()):
            return original
        prefix = original
        if prefix and not prefix.endswith("\n"):
            prefix += "\n"
        if prefix and not prefix.endswith("\n\n"):
            prefix += "\n"
        return prefix + "[sprkkr]\n" + "".join(
            f"{key} = {value}\n"
            for key, value in updates.items() if value is not None
        )

    written = set()
    replacement = []
    i = first + 1
    while i < last:
        line = lines[i]
        match = _OPTION.match(line.rstrip("\r\n"))
        key = match.group(2).strip().lower() if match else None
        if key not in updates:
            replacement.append(line)
            i += 1
            continue

        value = updates[key]
        if value is not None and key not in written:
            ending = "\r\n" if line.endswith("\r\n") else "\n"
            replacement.append(
                f"{match.group(1)}{match.group(2)}{match.group(3)}{value}{ending}"
            )
        written.add(key)
        i += 1
        # ConfigParser permits indented continuation lines. If the edited
        # option was multiline, the continuation belongs to that option.
        while (i < last and lines[i].strip() and lines[i][0].isspace()
               and not lines[i].lstrip().startswith(("#", ";"))):
            if _OPTION.match(lines[i].rstrip("\r\n")):
                break
            i += 1

    for key, value in updates.items():
        if key not in written and value is not None:
            if replacement and not replacement[-1].endswith(("\n", "\r")):
                replacement.append("\n")
            replacement.append(f"{key} = {value}\n")

    return "".join(lines[:first + 1] + replacement + lines[last:])


def write_ase_sprkkr(updates: dict[str, str | None], *, system: bool = False) -> None:
    """Atomically update ASE's [sprkkr] section and preserve other settings.

    Before changing anything, both the existing and proposed file must be
    valid INI. Preserve symlinks and permissions in the same way as the
    ASE2SPRKKR YAML writer.
    """
    path = ase_config_target(system=system).resolve()
    if not system and _is_system_config_path(path):
        raise PermissionError("System ASE configuration cannot be modified without --system")
    original = path.read_text(encoding="utf-8") if path.is_file() else ""
    parser = configparser.ConfigParser(
        interpolation=configparser.ExtendedInterpolation()
    )
    if original:
        parser.read_string(original, source=str(path))
    revised = _rewrite_section(original, updates)
    check = configparser.ConfigParser(
        interpolation=configparser.ExtendedInterpolation()
    )
    if revised:
        check.read_string(revised, source=str(path))

    if revised == original:
        return

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False,
        ) as stream:
            temporary = Path(stream.name)
            stream.write(revised)
            try:
                mode = stat.S_IMODE(path.stat().st_mode)
            except FileNotFoundError:
                pass
            else:
                os.chmod(temporary, mode)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
