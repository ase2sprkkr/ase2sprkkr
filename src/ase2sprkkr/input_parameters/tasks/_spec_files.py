"""Shared declarative native KKRSPEC file specifications.

Each callback returns a filename or None when the file is not used.
"""


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


def data_spc(parameters):
    if parameters.TASK.TASK() == "SPLEED":
        # Native SPLEED always opens the SPC (sometimes header-only).
        return "{DATASET_ADSI}_data.spc"
    return "{DATASET_ADSI}_data.spc" if _single_polarization(parameters) else None


def data_pol1_spc(parameters):
    return None if _single_polarization(parameters) else "{DATASET_ADSI}_data_pol1.spc"


def data_pol2_spc(parameters):
    return None if _single_polarization(parameters) else "{DATASET_ADSI}_data_pol2.spc"


SPEC_INPUT_FILES = (input_inp, struc_inp)
SPEC_OUTPUT_FILES = (
    "{DATASET_ADSI}_SPEC.out",
    data_spc,
    data_pol1_spc,
    data_pol2_spc,
)
