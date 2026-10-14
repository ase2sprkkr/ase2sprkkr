"""Helpers for task-specific auxiliary file manifests."""

import os

from ..task_definition import TaskDefinition


def resolve_file(directory, filename):
    """Resolve a task file name relative to the calculation directory."""
    filename = os.fspath(filename)
    return filename if os.path.isabs(filename) else os.path.join(directory, filename)


def _optional_value(section, name):
    try:
        return section[name]()
    except KeyError:
        return None


def tau_file(parameters, directory):
    """Return the effective native TAU-directory file path."""
    filename = _optional_value(parameters.TAU, "TAUDIR")
    if not filename:
        filename = f"{parameters.CONTROL.DATASET()}.tau.dir"

    filename = os.fspath(filename)
    if os.path.isabs(filename):
        return filename

    taupath = _optional_value(parameters.TAU, "TAUPATH")
    base = resolve_file(directory, taupath) if taupath else directory
    return os.path.join(base, filename)


class RequiredTauTaskDefinition(TaskDefinition):
    """Task requiring an existing TAU-directory file."""

    def input_files(self, parameters, directory):
        return (tau_file(parameters, directory), *super().input_files(parameters, directory))


class TauCacheTaskDefinition(TaskDefinition):
    """Task optionally reusing a TAU-directory cache."""

    def cache_files(self, parameters, directory):
        return (tau_file(parameters, directory), *super().cache_files(parameters, directory))
