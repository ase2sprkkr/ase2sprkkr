"""Option and section types used to expose task results."""

import os
import platform
import shutil
import subprocess
import functools
from typing import Any, Tuple

from ..common.container_definitions import SectionDefinition
from ..common.generated_configuration_definitions import GeneratedValueDefinition
from ..common.options import Option
from ..common.configuration_containers import Section
from ..common.file_utils import FilePath
from ..common.value_definitions import ValueDefinition
from ..common.warnings import ValidationResult
from ..result_files import ResultFile, read_plottable_file


class OutputOption(Option):
    @property
    def display_name(self):
        return getattr(self._definition, "display_name", self.name)

    def value_label(self):
        return self()

    def actions(self):
        return ("plot",) if callable(getattr(self, "plot", None)) else ()


class DataOutputOption(OutputOption):
    def value_label(self):
        return "<data...>"

    def actions(self):
        return ("data",)

    def data(self):
        return self()


class OutputFileOption(OutputOption):
    _FILE_DISPLAY_NAMES = {
        "output": "Output file",
        "input": "Input file",
        "potential": "Potential file",
        "converged": "Converged potential file",
        "sfn": "Shape functions (SFN)",
        "dos": "Density of states (DOS)",
        "bloch-sf": "Bloch spectral function (BSF)",
        "spc": "Photoemission spectrum",
        "jxc": "Exchange couplings",
        "dij": "Exchange tensors",
        "dmi": "Dzyaloshinskii–Moriya vectors",
    }

    _TYPE_DISPLAY_NAMES = {
        "dos": "Density of states (DOS)",
        "bsf": "Bloch spectral function (BSF)",
        "spc": "Photoemission spectrum",
        "jxc": "Exchange-interaction data",
        "rat": "X-ray spectroscopy data",
        "sfn": "Shape functions (SFN)",
        "agr": "Grace plot (AGR)",
    }

    _FILE_INFO = {
        "output": "Main text output produced by SPR-KKR.",
        "input": "Input parameters used for this calculation.",
        "potential": "Potential used as input for this calculation.",
        "converged": "Potential produced by the converged SCF calculation.",
        "sfn": "Shape-function data used by a full-potential calculation.",
        "agr": "Grace plot project containing one or more datasets.",
    }

    _TYPE_INFO = {
        "dos": "Density-of-states data produced by the calculation.",
        "bsf": "Bloch spectral-function data produced by the calculation.",
        "spc": "Photoemission spectrum data produced by the calculation.",
        "jxc": "Exchange-interaction data produced by the calculation.",
        "rat": "X-ray spectroscopy data produced by the calculation.",
        "sfn": "Shape-function data used by a full-potential calculation.",
    }

    @property
    def info(self):
        return (
            self._TYPE_INFO.get(self.file_type)
            or self._FILE_INFO.get(self.name.casefold())
            or f"File produced or used by the calculation: {self.name}."
        )

    @property
    def display_name(self):
        return (
            self._FILE_DISPLAY_NAMES.get(self.name.casefold())
            or self._TYPE_DISPLAY_NAMES.get(self.file_type)
            or super().display_name
        )

    def path(self):
        owner = getattr(self._container, "result", None)
        return owner.path_to(self.name) if owner is not None else os.fspath(self())

    @property
    def file_type(self):
        extension = os.path.splitext(self())[1].lstrip(".").lower()
        # Grace project format is authoritative even for a named output option
        # with an older task-specific file type hint.
        if extension == "agr":
            return "agr"
        if self._definition.file_type:
            return self._definition.file_type
        return extension if extension in ResultFile.definitions else None

    def definition(self):
        return ResultFile.definitions[self.file_type]

    def value_label(self):
        return f"<{self.file_type} file>" if self.file_type else str(self())

    def actions(self):
        # Parsed paths can describe missing optional or mandatory outputs.
        # Do not offer file operations before the actual file is accessible.
        if not os.path.isfile(self.path()):
            return ()
        out = ["open", "open_directory", "save"]
        if self.file_type == "agr":
            out.append("plot")
        elif self.file_type:
            result_class = self.definition().result_class
            if result_class.can_be_plotted():
                out.append("plot")
            for action in result_class.additional_actions:
                # Several unrelated formats use .spc. When the filename only
                # supplies that ambiguous extension, do not advertise ARPES
                # conversion for SPLEED/BAND/AIPES calculations.
                if action == "export_nexus" and self.file_type == "spc":
                    owner = getattr(self._container, "result", None)
                    if str(getattr(owner, "task_name", "")).upper() != "ARPES":
                        continue
                out.append(action)
        return tuple(out)

    def __getattr__(self, name):
        if not name.startswith("_") and self.file_type:
            result_class = self.definition().result_class
            if name in result_class.additional_actions:
                return functools.partial(self._run_file_action, name)
        raise AttributeError(f"{type(self).__name__!s} has no attribute {name!r}")

    def _run_file_action(self, action, parent=None, **kwargs):
        return getattr(self.parsed_file(), action)(**kwargs)

    def parsed_file(self, potential=None):
        """Parse the file, using its task result as context when available."""
        owner = self._container.result
        if owner is not None:
            return owner.parse_result_file(self, potential=potential)
        return read_plottable_file(
            self.path(),
            try_only=self.file_type,
            potential=potential,
        )

    def export_nexus(self, parent=None):
        """Export an eligible ARPES spectrum using the canonical NeXus writer."""
        import h5py  # noqa: F401 - check optional dependency before file dialog
        from PyQt6.QtWidgets import QFileDialog

        spectrum = self.parsed_file()
        exporter = getattr(spectrum, "to_nexus", None)
        if not callable(exporter):
            raise ValueError(f"NeXus export is unavailable for {self.path()}")

        proposed = os.path.splitext(os.fspath(self.path()))[0] + ".nxs"
        filename, _ = QFileDialog.getSaveFileName(
            parent, "Export ARPES Spectrum to NeXus", proposed,
            "NeXus files (*.nxs);;All files (*)",
        )
        if not filename:
            return None
        return exporter(filename)

    def plot(self, parent=None):
        self.parsed_file().plot()

    def save(self, parent=None):
        from PyQt6.QtWidgets import QFileDialog

        extension = (
            self.definition().extension
            if self.file_type and self.file_type != "agr"
            else os.path.splitext(self.path())[1].lstrip(".") or "*"
        )
        filename, _ = QFileDialog.getSaveFileName(parent, "Save File", "", f"(*.{extension});;All Files (*)")
        if filename:
            shutil.copy(self.path(), filename)

    def open(self, parent=None):
        if platform.system() == "Windows":
            os.startfile(self.path())
        elif platform.system() == "Darwin":
            subprocess.run(["open", self.path()])
        else:
            subprocess.run(["xdg-open", self.path()])

    def open_directory(self, parent=None):
        directory = os.path.dirname(os.path.abspath(self.path()))
        if platform.system() == "Windows":
            subprocess.run(["explorer", directory])
        elif platform.system() == "Darwin":
            subprocess.run(["open", directory])
        else:
            subprocess.run(["xdg-open", directory])


class ResultFileDefinition(ValueDefinition):
    result_class = OutputFileOption
    name_in_grammar = False

    def __init__(self, name, file_type=None):
        self.file_type = file_type
        super().__init__(name, str, is_optional=True)

    def convert_value(self, option: Option, value: Any, item: bool = False) -> FilePath:
        """Resolve ``value`` relative to the result owning ``option``.

        ``item`` is accepted for the common value-definition conversion API.
        """
        owner = getattr(option._container, "result", None)
        directory = owner.directory if owner is not None else None
        return FilePath(value, directory)

    def validate_converted(
        self, option: Option, value: Any, why: str = "set", item: bool = False
    ) -> Tuple[ValidationResult, ...]:
        """Accept converted ``value``; remaining arguments match the common API."""
        return ()


class OutputFilesSection(Section):
    def add_file(self, name, filename, file_type=None):
        if name in self:
            option = self[name]
            option.set(filename)
            if file_type is not None:
                option._definition.file_type = file_type
            return option
        definition = ResultFileDefinition(name, file_type)
        definition.added_to_container(self._definition)
        self._definition._members[name] = definition
        option = definition.create_object(self)
        option.set(filename)
        self._add(option)
        return option

    def set_file_type(self, name, file_type):
        """Set the parser type for an already registered output file."""
        self[name]._definition.file_type = file_type
        return self[name]


def create_files_section(result=None):
    definition = SectionDefinition("FILES", [], result_class=OutputFilesSection)
    out = definition.create_object()
    object.__setattr__(out, "result", result)
    return out


def generated_output(
    name,
    getter,
    *,
    result_class=OutputOption,
    plot=None,
    info=None,
    display_name=None,
):
    definition = GeneratedValueDefinition(
        name,
        getter,
        result_class=result_class,
        plot=plot,
        info=info,
    )
    definition.display_name = display_name or name
    return definition


def create_output_section(name, members):
    return SectionDefinition(name, members).create_object()
