"""The reader for JXC task, providing exchange couplings, exchange tensors, and Dzyaloshinski-Moriya vectors."""

import warnings
import re

from ..task_result import TaskResult, KkrOutputReader
from ..sprkkr_output_reader import SprKkrOutputParser
from ...common.decorators import cached_property
from ...bindings.uppasd import write_pos_file, write_mom_file, write_inpsd_file
from ...common.process_output_reader import readline_until


class JxcResult(TaskResult):
    """JXC task provides access to the generated JXC output file.
    It provides access to exchange couplings, exchange tensors, and Dzyaloshinski-Moriya vectors,
    using the :py:attr:`~jxc`, :py:attr:`~dij` and :py:attr:`~dmi` properties.
    Not all of them may be present in the output, depending on the settings of the JXC task.
    """

    @cached_property
    def dmi_filename(self):
        """Dzyaloshinski-Moriya vectors file name"""
        try:
            return self.path_to("dmi")
        except KeyError:
            raise AttributeError("Dzyaloshinski-Moriya vectors file not found in output")

    @cached_property
    def dmi(self):
        """Dzyaloshinski-Moriya vectors."""
        return self.output_values["dmi"].parsed_file()

    @cached_property
    def dij_filename(self):
        """Exchange tensors file name"""
        try:
            return self.path_to("dij")
        except KeyError:
            raise AttributeError("Exchange tensors file not found in output")

    @cached_property
    def dij(self):
        """Exchange tensors."""
        return self.output_values["dij"].parsed_file()

    @cached_property
    def jxc_filename(self):
        """Exchange couplings file name"""
        try:
            return self.path_to("jxc")
        except KeyError:
            raise AttributeError("Exchange couplings file not found in output")

    @cached_property
    def jxc(self):
        """Exchange couplings."""
        return self.output_values["jxc"].parsed_file()

    @property
    def output_values(self):
        files = {"jxc": "Exchange couplings", "dij": "Exchange tensors", "dmi": "Dzyaloshinski-Moriya vectors"}
        for key in files:
            if key in self.files:
                self.files.set_file_type(key, "jxc")
        return super().output_values

    def write_uppasd_files(
        self,
        directory=None,
        jxc_file="jfile.dat",
        dmi_file="dmfile.dat",
        pos_file="posfile.dat",
        mom_file="momfile.dat",
        inpsd_file="inpsd.dat",
        **kwargs,
    ):
        """Write the computed exchange couplings, exchange tensors and Dzyaloshinski-Moriya vectors
        to UppASD input files in the specified directory."""
        output_file = None
        args = {"dmi": "biqdm", "jxc": "exchange"}
        for key in ["jxc", "dmi"]:
            if key in self.files:
                file_name = locals()[key + "_file"]
                if file_name is False:
                    continue
                output_file = self.output_values[key].parsed_file()
                output_file.write_uppasd_file(file_name=file_name, directory=directory)
                argname = args[key]
                if argname in kwargs:
                    warnings.warn(
                        f"Option `{argname}` already present in the kwargs. I have no place to wrute the name of {key} file into inpsd.dat"
                    )
                else:
                    kwargs[argname] = file_name

        if pos_file is not False:
            write_pos_file(self.atoms, pos_file, jxc_ouput_file=output_file, directory=directory)
        if mom_file is not False:
            write_mom_file(self.atoms, mom_file, jxc_ouput_file=output_file, directory=directory)
        if inpsd_file is not False:
            write_inpsd_file(self.atoms, posfile=pos_file, momfile=mom_file, **kwargs)

    def run_uppasd(
        self,
        executable="uppasd",
        directory=None,
        jxc_file="jfile.dat",
        dmi_file="dmfile.dat",
        pos_file="posfile.dat",
        mom_file="momfile.dat",
        inpsd_file="inpsd.dat",
        **kwargs,
    ):
        self.write_uppasd_files(
            **{k: v for k, v in locals().items() if k not in ("self", "kwargs", "executable")}, **kwargs
        )
        import subprocess
        from ...common.file_utils import Directory

        with Directory(directory).chdir():
            subprocess.run([executable])


class JxcOutputParser(SprKkrOutputParser):
    async def read_output(self, stdout, result):
        await self.read_commons(stdout, result)

        pattern = re.compile(
            rb"\s*(?:(?:"
            + rb")|(?:".join(
                [
                    rb"(?P<file>results written to file:(?P<file_name>.+_(?P<file_type>XCPLTEN_Jij|XCPLTEN_Dij|DMIVEC_Dij)\.dat|.*))",
                    rb"Curie temperature within mean field approximation  T_C =\s*(?P<curie_temp>[0-9.]+) K",
                ]
            )
            + rb"))\s*"
        )
        file_types = {b"XCPLTEN_Jij": "jxc", b"XCPLTEN_Dij": "dij", b"DMIVEC_Dij": "dmi"}

        match = None

        def important(line):
            nonlocal match
            match = pattern.match(line)
            return match

        result.mean_field_curie_temperature = None
        while True:
            line = await readline_until(stdout, important)
            if line is None:
                break
            if match.group("file"):
                file_type = file_types.get(match.group("file_type"), None)
                if file_type is not None:
                    result.files.add_file(file_type, match.group("file_name").decode("utf8"), "jxc")
                else:
                    warnings.warn(f"Unexpected file type in JXC output: {line}")
            if match.group("curie_temp"):
                result.mean_field_curie_temperature = float(match.group("curie_temp"))

        if result.mean_field_curie_temperature is None:
            warnings.warn("Curie temperature not found in JXC output")


class JxcOutputReader(KkrOutputReader):
    result_class = JxcResult
    parser_class = JxcOutputParser
