"""Shared native KKRSPEC file specifications.

Callbacks return fully resolved filenames or None when unused.
"""


from ..task_definition import format_file_template


def _spc_file(parameters, name):
    return format_file_template(parameters, "{dataset}_" + name)


def input_inp(parameters):
    return "input.inp" if parameters.TASK.INPVER() == 0 else None


def struc_inp(parameters):
    return "struc.inp" if parameters.TASK.STRVER() == 0 else None


def _single_polarization(parameters):
    # BAND/SPLEED never enter the photon-polarization output branch.
    if parameters.TASK.TASK() in ("BAND", "SPLEED"):
        return True

    # NPOL selects the actual number of native outputs; POL_P specifies
    # the polarization orientation, not the number of calculated channels.
    # NPOL=0 (unpolarized/dichroism) and NPOL=3 (dichroism) use two.
    return str(parameters.SPEC_PH.NPOL()) not in ("0", "3")


def spleed_spc(parameters):
    """SPLEED always opens a structured SPC file in its default NPOL=1 mode.

    SPEC_RSLAB calls WRHEAD independently of TYP and NREL. For TYP=0 the
    file contains only a header; its two-column diagram is in SPEC.out.
    Modern SPEC_INPUT uses NREL=1, and TYP=1/2 write seven data columns.
    """
    return _spc_file(parameters, "data.spc") if parameters.TASK.TASK() == "SPLEED" else None


def data_spc(parameters):
    if parameters.TASK.TASK() == "SPLEED":
        return spleed_spc(parameters)
    return _spc_file(parameters, "data.spc") if _single_polarization(parameters) else None


def data_pol1_spc(parameters):
    return None if _single_polarization(parameters) else _spc_file(parameters, "data_pol1.spc")


def data_pol2_spc(parameters):
    return None if _single_polarization(parameters) else _spc_file(parameters, "data_pol2.spc")


SPEC_INPUT_FILES = (input_inp, struc_inp)
SPEC_OUTPUT_FILES = (
    "{dataset}_SPEC.out",
    data_spc,
    data_pol1_spc,
    data_pol2_spc,
)
