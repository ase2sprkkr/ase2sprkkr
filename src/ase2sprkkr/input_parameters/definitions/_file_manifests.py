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


class RequiredTauInputParametersDefinition(InputParametersDefinition):
    """Task definition requiring an existing TAU-directory file."""

    def input_files(self, parameters, directory):
        return (tau_file(parameters, directory), *super().input_files(parameters, directory))


class TauCacheInputParametersDefinition(InputParametersDefinition):
    """Task definition optionally reusing a TAU-directory cache."""

    def cache_files(self, parameters, directory):
        return (tau_file(parameters, directory), *super().cache_files(parameters, directory))


class NrSpectroscopyInputParametersDefinition(InputParametersDefinition):
    """Shared native output contract for NRAPS/NRAES convolution tasks."""

    def __init__(self, *args, output_extension, broadened_dos=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.output_extension = output_extension
        self.broadened_dos = broadened_dos

    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        dataset = parameters.CONTROL.DATASET()
        extension = self.output_extension

        # Native naming may insert ADSI / atom / shell information before the
        # extension, so retain that dynamic middle part as a constrained glob.
        files.append(resolve_file(directory, f"{dataset}*.{extension}"))

        if parameters.TASK.CALCLS():
            files.append(resolve_file(directory, f"{dataset}*.{extension}.??"))

        if parameters.TASK.OUTGNU():
            # L1,L2 are one-character angular-momentum labels; S1,S2 are U/D.
            files.extend(resolve_file(directory, f"{prefix}??[UD][UD]")
                         for prefix in ("D", "E", "DE", "W"))

        if self.broadened_dos:
            gamma = _optional_value(parameters.TASK, "GAMMA")
            if gamma is not None and gamma > 0:
                files.append(resolve_file(directory, f"{dataset}.dos.brd"))

        return tuple(files)
