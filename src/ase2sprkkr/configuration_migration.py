"""Migration helpers for the deprecated executable-Python user configuration."""

from __future__ import annotations

import ast
import importlib.machinery
import os
import pickle
import shutil
import subprocess
import sys
import tempfile
import types
from pathlib import Path


def load_legacy_configuration(file, definition):
    """Execute a legacy hook in a fresh process and import its explicit values.

    Both public config aliases refer to the child configuration. Arbitrary
    imports, registry changes and failed hooks cannot mutate the caller.
    """
    from .configuration_yaml import _set_persisted_configuration

    source = Path(file).resolve()
    with tempfile.TemporaryDirectory(prefix="ase2sprkkr-migrate-") as directory:
        schema = Path(directory) / "definition.pickle"
        snapshot = Path(directory) / "values.pickle"
        with schema.open("wb") as stream:
            pickle.dump(definition, stream)
        environment = dict(os.environ, ASE2SPRKKR_NO_USER_PROFILE="1")
        result = subprocess.run(
            [sys.executable, "-c",
             ("from ase2sprkkr.configuration_migration import _legacy_worker; "
              "import sys; _legacy_worker(*sys.argv[1:])"),
             os.fspath(source), os.fspath(schema), os.fspath(snapshot)],
            env=environment, stderr=subprocess.PIPE, text=True, check=False,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"Loading legacy configuration {source} failed in its worker "
                f"(exit {result.returncode}):\n{result.stderr.strip()}"
            )
        with snapshot.open("rb") as stream:
            values = pickle.load(stream)
    migrated = definition.create_object()
    migrated.profiles._legacy_write_warning_emitted = True
    _set_persisted_configuration(migrated, values)
    return migrated


def _legacy_worker(source, schema, snapshot):
    """Child-only entry point; exchange pickles only inside the private tempdir."""
    import ase2sprkkr

    from . import configuration

    with open(schema, "rb") as stream:
        migrated = pickle.load(stream).create_object()
    migrated.profiles._legacy_write_warning_emitted = True
    configuration.config = ase2sprkkr.config = migrated
    loader = importlib.machinery.SourceFileLoader("ase2sprkkr.personal_migration", source)
    module = types.ModuleType(loader.name)
    loader.exec_module(module)
    values = migrated.as_dict(only_changed="explicit", generated=False) or {}
    with open(snapshot, "wb") as stream:
        pickle.dump(values, stream)


def _migrated_profile_name(profiles):
    """Return a free conventional name for migrated execution settings."""
    name = "local"
    index = 2
    while name in profiles:
        name = f"local-{index}"
        index += 1
    return name


def _convert_legacy_execution_profile(config):
    """Move non-default deprecated execution settings to a normal profile."""
    profiles = config.profiles
    legacy = profiles._legacy_profile()
    if legacy is None or not profiles._legacy_available():
        return None

    active_profile = profiles.resolve_name()
    changed = legacy.as_dict(only_changed=True, generated=False) or {}
    values = {
        name: changed[name]
        for name in profiles.legacy_fields
        if name in changed
    }
    name = _migrated_profile_name(profiles)

    profiles._values.pop("legacy", None)
    profile = profiles.add(name)
    profile.set(values, unknown="fail")

    # Preserve profile precedence.  Without ASE the legacy profile was the
    # fallback and must become explicit.  If ASE is available it already won
    # over legacy, so leave default_profile unset and keep ASE active.
    if active_profile == "legacy":
        config.default_profile.set(name)
    return name


def migrate_legacy_configuration(source, destination, definition):
    """Convert a legacy Python config file to YAML and back up the source."""
    from .configuration_yaml import save_configuration_yaml

    source = Path(source)
    backup = source.with_name(source.name + ".bak")
    shutil.copy2(source, backup)
    migrated = load_legacy_configuration(source, definition)
    _convert_legacy_execution_profile(migrated)
    save_configuration_yaml(migrated, destination)
    return migrated


def clean_legacy_configuration(file):
    """Remove migrated top-level ``config.* = ...`` statements from a script.

    Only plain top-level assignments are removed. Conditional or otherwise
    dynamic uses of ``config`` are deliberately retained. A directly preceding
    contiguous comment block is removed together with an assignment. The
    ``from ase2sprkkr.configuration import config`` import is removed only when
    no remaining code refers to ``config``. Runs of blank lines left by removed
    blocks are collapsed to a single blank line.
    """
    path = Path(file)
    source = path.read_text(encoding="utf8")
    tree = ast.parse(source, filename=os.fspath(path))
    lines = source.splitlines(keepends=True)

    assignments = [
        node for node in tree.body
        if isinstance(node, (ast.Assign, ast.AnnAssign))
        and _is_config_assignment(node)
    ]
    remove = set()
    for node in assignments:
        start = node.lineno - 1
        end = node.end_lineno or node.lineno
        while start > 0:
            previous = lines[start - 1]
            if previous.lstrip().startswith("#"):
                start -= 1
                continue
            break
        remove.update(range(start, end))

    remaining_source = "".join(
        line for index, line in enumerate(lines) if index not in remove
    )
    remaining_tree = ast.parse(remaining_source, filename=os.fspath(path))
    if not any(
        isinstance(node, ast.Name) and node.id == "config"
        for node in ast.walk(remaining_tree)
        if not _is_config_import(node)
    ):
        remove.update(_config_import_lines(tree))

    cleaned = _collapse_blank_lines(
        line for index, line in enumerate(lines) if index not in remove
    )
    path.write_text(cleaned, encoding="utf8")
    return len(assignments)


def _collapse_blank_lines(lines):
    result = []
    previous_blank = False
    for line in lines:
        blank = not line.strip()
        if blank and previous_blank:
            continue
        result.append(line)
        previous_blank = blank
    return "".join(result)


def _is_config_assignment(node):
    if isinstance(node, ast.Assign):
        targets = node.targets
    else:
        targets = (node.target,)
    return bool(targets) and all(_config_attribute(target) for target in targets)


def _config_attribute(node):
    while isinstance(node, ast.Attribute):
        node = node.value
    return isinstance(node, ast.Name) and node.id == "config"


def _is_config_import(node):
    return (
        isinstance(node, ast.ImportFrom)
        and node.module == "ase2sprkkr.configuration"
        and any(alias.name == "config" for alias in node.names)
    )


def _config_import_lines(tree):
    lines = set()
    for node in tree.body:
        if _is_config_import(node):
            lines.update(range(node.lineno - 1, node.end_lineno or node.lineno))
    return lines
