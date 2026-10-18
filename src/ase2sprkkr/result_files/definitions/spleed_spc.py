"""Native KKRSPEC SPLEED spectra in seven numeric columns."""
from ..result_files_definitions import ResultFileDefinition, ResultFileValueDefinition as V, create_result_file_definition
from ..result_files import CommonResultFile
from ...common.grammar_types import NumpyArray, Prefixed

class SPLEEDSPCResultFile(CommonResultFile):
    """Retain the seven native SPLEED columns without guessing their meaning."""

class SPLEEDSPCDefinition(ResultFileDefinition):
    result_class = SPLEEDSPCResultFile
    @property
    def extension(self):
        return "spc"

definition = create_result_file_definition(
    "SPLEED",
    [
        V("NT", int),
        V("NP", int),
        V("COMMENT", Prefixed("#"), name_in_grammar=False),
        V("RAW_DATA", NumpyArray(written_shape=(-1, 7)), name_in_grammar=False),
    ],
    name="SPLEED_SPC",
    cls=SPLEEDSPCDefinition,
    info="KKRSPEC SPLEED output: seven native numerical columns.",
)
