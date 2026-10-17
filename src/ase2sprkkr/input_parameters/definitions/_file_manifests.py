"""Helpers for task-specific auxiliary file manifests."""

import os

from ..input_parameters_definitions import InputParametersDefinition


def resolve_file(directory, filename):
    """Resolve a task file name relative to the calculation directory."""
    filename = os.fspath(filename)
    return filename if os.path.isabs(filename) else os.path.join(directory, filename)


def _optional_value(section, name):
    try:
        return section[name]()
    except KeyError:
        return None


def native_dataset(parameters):
    """Return SPR-KKR's effective ``DATSET`` value.

    Native executables preserve ``CONTROL.DATASET`` separately as ``DATSET0``
    and append ``_' + CONTROL.ADSI`` to ``DATSET`` when ADSI is non-empty.
    Most calculation-generated files use the latter value.
    """
    dataset = os.fspath(parameters.CONTROL.DATASET())
    adsi = _optional_value(parameters.CONTROL, "ADSI")
    return f"{dataset}_{adsi}" if adsi else dataset


def native_file(parameters, directory, suffix):
    """Resolve a file whose name is ``DATSET + suffix`` in SPR-KKR."""
    return resolve_file(directory, f"{native_dataset(parameters)}{suffix}")


def tau_file(parameters, directory):
    """Return the effective native TAU-directory path."""
    filename = _optional_value(parameters.TAU, "TAUDIR")
    if not filename:
        filename = f"{native_dataset(parameters)}.tau.dir"

    filename = os.fspath(filename)
    if os.path.isabs(filename):
        return filename

    taupath = _optional_value(parameters.TAU, "TAUPATH")
    base = resolve_file(directory, taupath) if taupath else directory
    return os.path.join(base, filename)


class RequiredTauInputParametersDefinition(InputParametersDefinition):
    """Task definition requiring an existing TAU-directory file."""

    def input_files(self, parameters, directory):
        return (tau_file(parameters, directory), *super().input_files(parameters, directory))


class TauCacheInputParametersDefinition(InputParametersDefinition):
    """Task definition optionally reusing a TAU-directory cache."""

    def cache_files(self, parameters, directory):
        return (tau_file(parameters, directory), *super().cache_files(parameters, directory))
