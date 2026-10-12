"""Migration helpers for the deprecated executable-Python user configuration."""

from __future__ import annotations

import ast
import importlib.machinery
import os
import shutil
import types
from pathlib import Path


def load_legacy_configuration(file, definition):
    """Execute one legacy config file against a fresh configuration object."""
    from . import configuration

    migrated = definition.create_object()
    original = configuration.config
    configuration.config = migrated
    try:
        loader = importlib.machinery.SourceFileLoader(
            "ase2sprkkr.personal_migration", os.fspath(file)
        )
        module = types.ModuleType(loader.name)
        loader.exec_module(module)
    finally:
        configuration.config = original
    return migrated


def migrate_legacy_configuration(source, destination, definition):
    """Convert a legacy Python config file to YAML and back up the source."""
    from .configuration_yaml import save_configuration_yaml

    source = Path(source)
    backup = source.with_name(source.name + ".bak")
    shutil.copy2(source, backup)
    migrated = load_legacy_configuration(source, definition)
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
