"""Option and section types used to expose task results."""

import os
import platform
import shutil
import subprocess

from ..common.container_definitions import SectionDefinition
from ..common.generated_configuration_definitions import GeneratedValueDefinition
from ..common.options import Option
from ..common.configuration_containers import Section
from ..common.value_definitions import ValueDefinition
from ..output_files.output_files import OutputFile


class OutputOption(Option):
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
    def path(self):
        owner = getattr(self._container, "result", None)
        return owner.path_to(self.name) if owner is not None else self()

    @property
    def file_type(self):
        return self._definition.file_type

    def definition(self):
        return OutputFile.definitions[self.file_type]

    def value_label(self):
        return f"<{self.file_type} file>" if self.file_type else self()

    def actions(self):
        out = ["open", "open_directory", "save"]
        if self.file_type and self.definition().result_class.can_be_plotted():
            out.append("plot")
        return tuple(out)

    def parsed_file(self):
        return OutputFile.from_file(self.path(), try_only=self.file_type)

    def plot(self, parent=None):
        self.parsed_file().plot()

    def save(self, parent=None):
        from PyQt6.QtWidgets import QFileDialog

        extension = self.definition().extension if self.file_type else "*"
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


class OutputFileDefinition(ValueDefinition):
    result_class = OutputFileOption
    name_in_grammar = False

    def __init__(self, name, file_type=None):
        self.file_type = file_type
        super().__init__(name, str, is_optional=True)


class OutputFilesSection(Section):
    def add_file(self, name, filename, file_type=None):
        if name in self:
            option = self[name]
            option.set(filename)
            if file_type is not None:
                option._definition.file_type = file_type
            return option
        definition = OutputFileDefinition(name, file_type)
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


def generated_output(name, getter, *, result_class=OutputOption, plot=None):
    return GeneratedValueDefinition(name, getter, result_class=result_class, plot=plot)


def create_output_section(name, members):
    return SectionDefinition(name, members).create_object()
