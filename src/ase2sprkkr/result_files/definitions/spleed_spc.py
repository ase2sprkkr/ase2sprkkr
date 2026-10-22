"""Structured seven-column relativistic SPLEED data from KKRSPEC.

The columns are taken from SPLOUT in SPR-KKR 9.7 (spec_rslabsub.f):
incident energy, angle, diffraction intensity, normalized spin asymmetry,
and three components of the outgoing spin polarization.

SPLOUT's two-column diagnostics go to NOUT1 (stdout), not this SPC file.
"""

from ..result_files_definitions import (
    ResultFileDefinition,
    ResultFileValueDefinition as V,
    create_result_file_definition,
)
from ..result_files import CommonResultFile
from ...common.grammar_types import NumpyArray, Prefixed
from ...common.generated_configuration_definitions import NumpyViewDefinition as NV


class SPLEEDSPCResultFile(CommonResultFile):
    """Relativistic SPLEED diffraction intensities and spin polarization."""


class SPLEEDSPCDefinition(ResultFileDefinition):
    result_class = SPLEEDSPCResultFile

    @property
    def extension(self):
        return "spc"


def column(index):
    return slice(None), index


definition = create_result_file_definition(
    "SPLEED",
    [
        V("NT", int),
        V("NP", int),
        V("COMMENT", Prefixed("#"), name_in_grammar=False),
        V("RAW_DATA", NumpyArray(shape=(-1, 7), written_shape=(-1, 7)), name_in_grammar=False),
        NV("ENERGY", "RAW_DATA", column(0), info="Incident electron kinetic energy (eV)"),
        NV("THETA", "RAW_DATA", column(1), info="Electron angle (degrees)"),
        NV("INTENSITY", "RAW_DATA", column(2), info="Diffracted beam intensity (RINTG)"),
        NV(
            "ASYMMETRY", "RAW_DATA", column(3),
            info="Spin asymmetry (I+ - I-)/(I+ + I-) divided by |P|",
        ),
        NV("P_X", "RAW_DATA", column(4), info="Outgoing spin polarization x"),
        NV("P_Y", "RAW_DATA", column(5), info="Outgoing spin polarization y"),
        NV("P_Z", "RAW_DATA", column(6), info="Outgoing spin polarization z"),
    ],
    name="SPLEED_SPC",
    cls=SPLEEDSPCDefinition,
    info="Seven-column relativistic KKRSPEC SPLEED diffraction output (SPLOUT).",
)
