"""Containers for configuration parameters of SPR-KKR task.

These configuration parameters are supplied to SPR-KKR executables as input files.
"""

from __future__ import annotations

import os
import io
import pkgutil
import importlib
import pyparsing as pp
from ..outputs.task_result import KkrOutputReader
from . import definitions
from ..sprkkr.configuration import ConfigurationFile, ConfigurationSection
from ..common.decorators import cached_class_property
from ..common.warnings import DataValidityError
from typing import Union
from ..configuration import config
from ..potentials.potentials import Potential
from ..common.doc import process_input_parameters_definition
from ..sprkkr.profile import SPRKKRProfile, executor_from_profile
from ..execution.persistence import resolve_execution_manager
import warnings


class InputSection(ConfigurationSection):
    """Input parameters sections has nothing special, yet."""


class InputParameters(ConfigurationFile):
    """It holds the configuration values for the task.

    According to its definition it can also prepare and execute the
    corresponding SPR-KKR calculation through an execution backend.
    """

    @property
    def potential(self):
        pot = getattr(self, "_potential", None)
        potfil = self.CONTROL.POTFIL.result
        if pot:
            if potfil == pot[0]:
                return pot[1]
        out = Potential.from_file(potfil) if potfil else None
        self._potential = potfil, out
        return out

    def resolve_executable_suffix(self, suffix: Union[str, bool]):
        """ " Return the suffix, that is appended after the name of SPR-KKR executable.

        Parameters
        ----------
        suffix
          - If str is given, it is left as is.
          - If True, return the default value:
                      config.executable.suffix
                     (content of SPRKKR_EXECUTABLE_SUFFIX env variable
                      or a user_defined_value)
          - If False, return ''
        """
        if suffix is False:
            return ""
        if suffix is True:
            return config.executables.suffix()
        return suffix

    def __init__(self, definition, inputfile=None, outputfile=False):
        super().__init__(definition)
        self._inputfile = inputfile
        self._outputfile = inputfile if outputfile is False else outputfile

    def copy(self, copy_values: bool = False):
        """Copy the parameters, preserving the originally requested task name."""
        out = super().copy(copy_values=copy_values)
        if hasattr(self, "_requested_task_name"):
            out._requested_task_name = self._requested_task_name
        return out

    @property
    def task_name(self):
        """Return the task name, as defined in the definition of the parameters
        (see InputParametersDefinition class)"""
        return self._definition.name

    @property
    def _execution_task_name(self):
        """Return the task encoded by the concrete input parameters."""
        task_section = getattr(self, "TASK", None)
        task_option = getattr(task_section, "TASK", None) if task_section is not None else None
        return task_option() if task_option is not None else self.task_name

    def uses_separate_single_site_contour(self):
        """Whether SPR-KKR evaluates single-site terms on a separate contour.

        The effective mode can be requested explicitly in ENERGY or CONTROL,
        is implied by the relativistic hyperfine-field decomposition, and is
        intrinsic to COMPTON calculations.  Missing task-specific options are
        treated as disabled.
        """
        switches = (
            ("ENERGY", "SPLITSS"),
            ("CONTROL", "SPLITSS"),
            ("CONTROL", "FSOHFF"),
        )
        for section, name in switches:
            try:
                if self[section][name]():
                    return True
            except (KeyError, AttributeError):
                pass
        return self.task_name.upper() == "COMPTON"

    def is_mpi(self, mpi=True):
        """Will be this task (the task described by this input parameters)
        runned using mpi?

        Parameters
        ----------
        mpi: list or str or bool
          Optional parameter. If False is given, return False regardless the type of the
          task, otherwise it is ignored. See InputParameters.mpi_runner for the explanation of the
          parameter.

        Return
        ------
        is_mpi: bool
          Will be this task runned using mpi?
        """
        return mpi and self._definition.mpi

    def input_files(self, directory):
        """Return task-specific input files for this concrete calculation."""
        return self._definition.input_files(self, directory)

    def output_files(self, directory):
        """Return task-specific output files for this concrete calculation."""
        return self._definition.output_files(self, directory)

    def cache_files(self, directory):
        """Return task-specific reusable cache files for this calculation."""
        return self._definition.cache_files(self, directory)

    def run_process(
        self,
        calculator,
        input_file,
        output_file,
        directory=".",
        print_output=None,
        read_callback=None,
        run_async=False,
        callback=None,
        gdb=False,
        profile=None,
        executor=None,
        persistent=None,
    ):
        """Run the SPR-KKR task described by these input parameters.

        ``run_async=False`` executes synchronously and returns a TaskResult.
        ``run_async=True`` returns an unstarted Runner.

        ``persistent`` may be ``False``/``None`` (disabled), ``True`` (use the
        process-wide default execution manager), or an explicit
        :class:`~ase2sprkkr.execution.ExecutionManager`.
        """
        print_output = calculator.value_or_default("print_output", print_output)
        reader = self.output_reader(calculator, print_output, directory)
        explicit_profile = profile is not None
        profile = profile or getattr(calculator, "profile", None) or SPRKKRProfile()
        if executor is None and not explicit_profile:
            executor = getattr(calculator, "executor", None)
        if executor is None:
            executor = executor_from_profile(profile)
        if persistent is None:
            persistent = getattr(calculator, "persistent", False)
        manager = resolve_execution_manager(persistent)

        execution_profile = profile
        if self._definition.mpi:
            execution_profile = profile.with_overrides(
                mpi=executor.resolve_mpi(profile.mpi)
            )

        try:

            def cleanup(out):
                input_file.close()
                if callback:
                    callback(out)

            potential_file = self.CONTROL.POTFIL.result
            if potential_file and not os.path.isabs(potential_file):
                potential_file = os.path.join(directory, potential_file)

            input_files = tuple(
                file
                for file in (
                    input_file,
                    potential_file,
                    *self.input_files(directory),
                )
                if file is not None
            )
            output_files = tuple(
                file
                for file in (
                    output_file,
                    *self.output_files(directory),
                )
                if file is not None
            )
            cache_files = tuple(
                file for file in self.cache_files(directory) if file is not None
            )

            spec = execution_profile.create_execution(
                self._definition.executable,
                input_file,
                output_file,
                directory,
                input_files=input_files,
                output_files=output_files,
                cache_files=cache_files,
                mpi_supported=self._definition.mpi,
                gdb=gdb,
            )
            if manager is None:
                runner = executor.create_run(reader, spec, callback=cleanup)
            else:
                runner = manager.create_run(
                    executor,
                    reader,
                    spec,
                    task=self._execution_task_name,
                    dataset=self.CONTROL.DATASET.result,
                    callback=cleanup,
                )
            if read_callback:
                runner.add_output_listener(read_callback)
            if run_async:
                return runner
            return runner.run()
        except FileNotFoundError as e:
            add = (
                "Cannot find SPR-KKR executable. Check executable_dir and executable_suffix "
                "in the selected SPR-KKR profile.\n"
            )
            if execution_profile.mpi:
                add += (
                    "Also, please check the MPI setting in the selected profile, or explicitly "
                    "disable MPI with mpi=False in SPRKKR.calculate().\n\n"
                )
            e.strerror = add + "SPRKKR cannot be run due to the following error: \n" + e.strerror
            raise

    def output_reader(self, calculator=None, print_output=False, directory=None):
        """Return the task-specific output reader/parser."""
        cls = self._definition.result_reader

        if not directory and calculator and calculator.directory:
            directory = calculator.directory

        if cls is None:
            task = self.TASK.TASK().lower()
            cls = KkrOutputReader.class_for_task(task)
        return cls(self, calculator, directory, print_output)

    def read_output_from_file(self, filename, directory=None):
        """Read and parse output of a previously run task."""
        directory = directory or os.path.dirname(filename)
        return self.output_reader(directory=directory).read_from_file(filename)

    def executable_params(self, directory=None):
        """
        Return
        ------
        ([executable], stdin, stdout) params for process.Popen
        """

    def __str__(self):
        return f"<Configuration container {self._get_path()}>"

    def set_option(self, name, value):
        for i in self._members:
            if name in i:
                self._members[i][name] = value
                break
        else:
            raise KeyError("No option with name {} in any of the members".format(name))

    @cached_class_property
    def definition_modules():
        names = (i for i in pkgutil.iter_modules(definitions.__path__))
        im = importlib.import_module
        modules = (im(".definitions." + i.name, __package__) for i in names)
        return {m.__name__.rsplit(".", 1)[-1].upper(): m for m in modules if hasattr(m, "input_parameters")}

    _definitions = {}

    @classmethod
    def definition(cls, name):
        name = name.upper()
        if name not in cls._definitions:
            if name == 'BSFEK' or name == 'BSFKK':
                return cls.definition('bsf')
            module = cls.definition_modules[name]
            ip = module.input_parameters
            if isinstance(ip, ipdefs.InputParametersDefinition):
                warnings.warn(
                    f"In the module '{module.__name__}' in file '{module.__file__}' "
                    "is an InputParametersDefinition object directly defined. "
                    "Consider to enclose its definition in a lambda function to "
                    "speed up ase2sprkkr loading."
                )
            else:
                try:
                    ip = ip()
                except Exception as e:
                    raise RuntimeError(
                        f"Can not load file {module.__name__} in {module.__file__} due to the following error:\n {e}"
                    ) from e
            cls._definitions[name] = ip
            process_input_parameters_definition(module, ip)
        return cls._definitions[name]

    @cached_class_property
    def definitions():
        # user = os.path.join(platformdirs.user_config_dir('ase2sprkkr', 'ase2sprkkr'), 'input_parameters')
        ip = InputParameters
        return {n: ip.definition(n) for n in ip.definition_modules}

    @classmethod
    def is_it_an_input_parameters_name(cls, name):
        name = name.upper()
        return name if name in cls.definition_modules else False

    @classmethod
    def create_input_parameters(cls, arg):
        """
        Create input_parameters object

        Parameters
        ----------
        arg: str or InputParameters
          If an InputParameters object is given, it is returned as is.
          If a string is given, it is interpreted either as a filename
          (from which the parameters are read) or the task name, for
          which the default parameters are used

        Return
        ------
        input_parameters: InputParameters
        """
        if isinstance(arg, str):
            name = cls.is_it_an_input_parameters_name(arg)
            if name:
                return cls.create(arg)
            return cls.from_file(arg)
        if isinstance(arg, io.IOBase):
            return cls.from_file(arg)
        return arg

    @classmethod
    def create(cls, name):
        """Create input parameters for the given task name
        Parameters
        ----------
        name: str
          Name of the task (e.g. 'SCF', 'PHAGEN')

        Return
        ------
        input_parameters: InputParameters
          Input parameters with the default values for the task.
        """
        out = cls(cls.definition(name))
        out._requested_task_name = name.upper()
        return out

    @classmethod
    def default_parameters(cls):
        """Create default input parameters"""
        return cls.create("SCF")

    @classmethod
    def from_file(cls, filename, allow_dangerous: bool = False):
        """Read an input file and create a new InputParameters object from the readed stuff

        Parameters
        ----------
        filename: str or file
          Input file (either its filename, or an open file)

        allow_dangerous
          Allow to read dangerous values of options: i.e. the values that do not fullfil the
          required type of the given option or its other requirements.
        """
        definitions = cls.definitions
        for d in definitions.values():
            try:
                out = d.read_from_file(filename, allow_dangerous=allow_dangerous)
                return out
            except (pp.ParseBaseException, DataValidityError) as e:
                last = e
        raise last

    def calculate(self, *args, **kwargs):
        """Create a calculator and run the input_parameters. See SPRKKR.calculate for the arguments"""
        calc = calculator.SPRKKR()
        calc.calculate(input_parameters=self, *args, **kwargs)

    def __repr__(self):
        d = self._definition
        out = d.configuration_type_name
        out = out + " for task " + d.name.upper()
        return out

    def change_task(self, task, retain_values=False):
        """Change the task to the given task. Retain the value of the options,
        that are present in the new task.
        """
        if retain_values:
            vals = self.to_dict(only_changed=True)
        self._definition = self.definition(task)
        self._requested_task_name = task.upper()
        self._init_members_from_the_definition()
        if retain_values:
            self.set(
                vals,
                unknown="ignore",
                retain_invalid="none",
                report_invalid="ignore",
            )

    def save_to_file(self, file, atoms=None, *, validate="save"):
        if self._definition.save_hook:
            self._definition.save_hook(getattr(file, "name", None), atoms)
        super().save_to_file(file, atoms, validate=validate)


# at least, to avoid a circular import
from ..sprkkr import calculator  # NOQA: E402
from . import input_parameters_definitions as ipdefs  # NOQA: E402
