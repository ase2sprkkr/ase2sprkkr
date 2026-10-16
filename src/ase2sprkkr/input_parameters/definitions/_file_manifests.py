"""Helpers for task-owned auxiliary file manifests."""

import os


def resolve_file(directory, filename):
    """Resolve a file name relative to the calculation directory."""
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


def bsf_file(parameters, directory):
    """Return the native BSF result path for the concrete DATASET/ADSI."""
    dataset = parameters.CONTROL.DATASET()
    adsi = parameters.CONTROL.ADSI()
    if not dataset or not adsi:
        return None
    return resolve_file(directory, f"{dataset}_{adsi}.bsf")


def kkrspec_uses_legacy_input(parameters):
    """Whether KKRSPEC reads externally supplied ``input.inp`` files."""
    return parameters.TASK.INPVER() == 0


def kkrspec_uses_legacy_structure(parameters):
    """Whether KKRSPEC reads an externally supplied structure file."""
    return parameters.TASK.STRVER() == 0


def kkrspec_single_spc_output(parameters):
    """Whether KKRSPEC writes one unpolarized/single SPC result file."""
    return parameters.TASK.TASK() == "BAND" or parameters.SPEC_PH.NPOL() == 1


def kkrspec_polarized_spc_output(parameters):
    """Whether KKRSPEC writes the pair of polarization-resolved SPC files."""
    return not kkrspec_single_spc_output(parameters)


def kkrspec_task_file_templates():
    """Return the shared file declarations of ARPES-family KKRSPEC tasks."""
    from ..task_definition import FileSpec

    return (
        (
            FileSpec("input.inp", condition=kkrspec_uses_legacy_input),
            FileSpec("in_bar.inp", condition=kkrspec_uses_legacy_input),
            FileSpec("in_structur.inp", condition=kkrspec_uses_legacy_structure),
        ),
        (
            FileSpec("{dataset}_data.spc", condition=kkrspec_single_spc_output),
            FileSpec(
                "{dataset}_data_pol1.spc",
                condition=kkrspec_polarized_spc_output,
            ),
            FileSpec(
                "{dataset}_data_pol2.spc",
                condition=kkrspec_polarized_spc_output,
            ),
        ),
    )
