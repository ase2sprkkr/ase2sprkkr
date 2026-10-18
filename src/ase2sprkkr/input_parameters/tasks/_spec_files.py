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

    photons = parameters.SPEC_PH
    polarization = photons.POL_P()
    if polarization in ("P", "C+"):
        return True
    if polarization in ("S", "C-"):
        return False
    return str(photons.NPOL()) == "1"


def data_spc(parameters):
    return "{dataset}_data.spc" if _single_polarization(parameters) else None


def data_pol1_spc(parameters):
    return None if _single_polarization(parameters) else "{dataset}_data_pol1.spc"


def data_pol2_spc(parameters):
    return None if _single_polarization(parameters) else "{dataset}_data_pol2.spc"


SPEC_INPUT_FILES = (input_inp, struc_inp)
SPEC_OUTPUT_FILES = (
    "{dataset}_SPEC.out",
    data_spc,
    data_pol1_spc,
    data_pol2_spc,
)
