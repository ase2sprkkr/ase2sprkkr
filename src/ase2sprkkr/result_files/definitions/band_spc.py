"""Native KKRSPEC BAND spectra in three numeric columns."""
from ..result_files_definitions import ResultFileDefinition, ResultFileValueDefinition as V, create_result_file_definition
from ..result_files import CommonResultFile
from ...common.grammar_types import NumpyArray, Prefixed
from ...common.generated_configuration_definitions import NumpyViewDefinition as NV

class BandSPCResultFile(CommonResultFile):
    """Raw band spectra. K_EIG's physical axis depends on NT."""

class BandSPCDefinition(ResultFileDefinition):
    result_class = BandSPCResultFile
    @property
    def extension(self):
        return "spc"

definition = create_result_file_definition(
    "BAND",
    [
        V("NT", int),
        V("NP", int),
        V("COMMENT", Prefixed("#"), name_in_grammar=False),
        V("RAW_DATA", NumpyArray(written_shape=(-1, 3)), name_in_grammar=False),
        NV("K_EIG", "RAW_DATA", (slice(None), 0)),
        NV("E_EIG", "RAW_DATA", (slice(None), 1)),
    ],
    name="BAND_SPC",
    cls=BandSPCDefinition,
    info="KKRSPEC BAND output: k coordinate, energy relative to EF, and auxiliary third column.",
)
