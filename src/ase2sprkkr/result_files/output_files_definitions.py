"""Internal compatibility imports for result-file format definitions."""

from .result_files_definitions import (
    BlankSeparator,
    OutputFileDefinition,
    OutputFileSection,
    OutputFileSectionDefinition,
    OutputFileValueDefinition,
    Separator,
    create_output_file_definition,
    output_file_header,
)

__all__ = [
    "BlankSeparator",
    "OutputFileDefinition",
    "OutputFileSection",
    "OutputFileSectionDefinition",
    "OutputFileValueDefinition",
    "Separator",
    "create_output_file_definition",
    "output_file_header",
]
