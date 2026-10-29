"""Helpers for task-specific auxiliary file manifests."""

import os

from ..input_parameters_definitions import InputParametersDefinition
from ._filenames import dataset_adsi


_XRAY_SHELLS = "KLMNO"
_XRAY_SUBSHELL_SUFFIXES = {"S": "1", "P": "23", "D": "45", "F": "67"}


def resolve_file(directory, filename):
    """Resolve a task file name relative to the calculation directory."""
    filename = os.fspath(filename)
    return filename if os.path.isabs(filename) else os.path.join(directory, filename)


def _optional_value(section, name):
    try:
        return section[name]()
    except KeyError:
        return None


def orbital_label(angular_momentum):
    """Return SPR-KKR's one-character FUNTXTL label for orbital angular momentum."""
    angular_momentum = int(angular_momentum)
    if angular_momentum <= 2:
        return "spd"[angular_momentum]
    return chr(ord("f") + angular_momentum - 3)


def xray_shell_label(core_level):
    """Return the native K/L/M/N/O filename label for a core shell."""
    principal = int(core_level[:-1])
    orbital = core_level[-1].upper()
    shell = _XRAY_SHELLS[principal - 1]
    if principal == 1:
        return shell
    return shell + _XRAY_SUBSHELL_SUFFIXES[orbital]


def core_level_output_name(parameters, extension, *, type_separator, shell_separator):
    """Return a constrained native filename pattern for one selected core shell."""
    dataset = dataset_adsi(parameters)
    shell = xray_shell_label(parameters.TASK.CL())
    return f"{dataset}{type_separator}*{shell_separator}{shell}.{extension}"


def tau_file(parameters, directory):
    """Return the effective native TAU-directory file path."""
    filename = _optional_value(parameters.TAU, "TAUDIR")
    if not filename:
        filename = dataset_adsi(parameters) + ".tau.dir"

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


class XrayRateInputParametersDefinition(TauCacheInputParametersDefinition):
    """XRAYSPEC definition with its native atom/core-shell rate filename."""

    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        name = core_level_output_name(
            parameters, "rat", type_separator="_", shell_separator="_"
        )
        files.append(resolve_file(directory, name))
        return tuple(files)


class NrSpectroscopyInputParametersDefinition(InputParametersDefinition):
    """Shared native output contract for NRAPS/NRAES convolution tasks."""

    def __init__(self, *args, output_extension, broadened_dos=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.output_extension = output_extension
        self.broadened_dos = broadened_dos

    def output_files(self, parameters, directory):
        files = list(super().output_files(parameters, directory))
        name = core_level_output_name(
            parameters, self.output_extension, type_separator=".", shell_separator="."
        )
        files.append(resolve_file(directory, name))

        if parameters.TASK.CALCLS():
            files.append(resolve_file(directory, f"{name}.??"))

        if parameters.TASK.OUTGNU():
            # D/E/W carry L1,L2 and the two-spin label; DE carries only L1,L2.
            files.extend(resolve_file(directory, f"{prefix}??[UD][UD]")
                         for prefix in ("D", "E", "W"))
            files.append(resolve_file(directory, "DE??"))

        if self.broadened_dos:
            gamma = _optional_value(parameters.TASK, "GAMMA")
            if gamma is not None and gamma > 0:
                files.append(resolve_file(directory, f"{dataset_adsi(parameters)}.dos.brd"))

        return tuple(files)
