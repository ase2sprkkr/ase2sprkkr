"""This module contains configuration, that could be changed, preferrably
by .config/ase2sprkkr/__init__.py file"""

import os
import re
from typing import Union
from .common.decorators import cache
from .common.grammar_types import CustomMixed, AlwaysQString, Array, Bool, Keyword, Integer
from .common.container_definitions import SectionDefinition, ConfigurationRootDefinition
from .common.configuration_containers import RootConfigurationContainer
from .common.value_definitions import ValueDefinition
from .common.generated_configuration_definitions import GeneratedValueDefinition
from .common.options import Option
from .profiles import ProfilesConfigurationContainer
import warnings
import shutil
import platformdirs


class Section(SectionDefinition):
    info_in_data_description = True
    dir_common_attributes = False


class ProfilesSection(Section):
    """Definition of the dict-like named execution-profile section."""

    repeated_class = ProfilesConfigurationContainer


def _get_suffix(*_):
    return os.environ.get("SPRKKR_EXECUTABLE_SUFFIX", "")


def _legacy_value(section, name, default=None):
    return section._container.profiles.legacy_value(name, default)


def _set_legacy_value(section, value, name):
    section._container.profiles.set_legacy_value(name, value)


def _legacy_mpi(section):
    return _legacy_value(section, "mpi")


def _set_legacy_mpi(section, value):
    _set_legacy_value(section, value, "mpi")


def _legacy_executable_dir(section):
    return _legacy_value(section, "executable_dir")


def _set_legacy_executable_dir(section, value):
    _set_legacy_value(section, value, "executable_dir")


def _legacy_executable_suffix(section):
    return _legacy_value(section, "executable_suffix", _get_suffix)


def _set_legacy_executable_suffix(section, value):
    _set_legacy_value(section, value, "executable_suffix")


def user_preferences_file():
    """Return filename with user preferences"""
    return os.path.join(platformdirs.user_config_dir("ase2sprkkr", "ase2sprkkr"), "__init__.py")


def load_user_preferences():
    """Load user defined preferences from
    ``$HOME/.config/ase2sprkkr/__init__.py``
    """
    file = user_preferences_file()

    try:
        if os.path.isfile(file):
            import types
            import importlib.machinery

            loader = importlib.machinery.SourceFileLoader("ase2sprkkr.personal", file)
            mod = types.ModuleType(loader.name)
            loader.exec_module(mod)
    except Exception as e:
        import warnings
        import traceback
        tb = "\n".join(traceback.format_exception(e))

        warnings.warn(
            f"Can not import {file} file with the user preferences:\n"
            f"{e}\n\n{tb}"
        )


@cache
def find_default_mpi_runner():
    if os.environ.get("SLURM_NODELIST") and shutil.which("srun"):
        return ["srun"]
    for r in ["mpirun", "mpirun.opmpirun", "mpirun.mpich", "mpiexec"]:
        if shutil.which(r):
            return [r]
    return False


@cache
def get_default_mpi_runner():

    out = find_default_mpi_runner()
    if out:
        return out
    if config.mpi_warning():
        warnings.warn("No MPI runner found. Disabling MPI!!!")


def mpi_runner(mpi):
    """Return a shell command to execute a mpi task.

    Parameters
    ----------
    mpi_runner: Union[bool,str,list,int]


      - If True is given, return the default mpi-runner
      - If False is given, no mpi-runner is returned.
      - If 'auto' is given, it is the same as True, however
           * no warning is given if no mpi is found
           * MPI is not used, if only one CPU is available
      - If a string is given, it is interpreted as a list of one item
      - If a list (of strings) is given, the user specified its own runner, use it as is
        as the parameters for subprocess.run.
      - If an integer is given, it is interpreted as the number of
        processes: the default mpi-runner is used, and the parameters
        to specify the number of processes.

    Return
    ------
    mpi_runner: list
      List of strings with the executable and its parameters, e.g.

      ::

          ["mpirun", "-np", "4"]
    """
    if mpi is None:
        mpi = config.running.mpi()
    if mpi is False:
        return None
    if mpi is True:
        mpi = find_default_mpi_runner()
    if isinstance(mpi, list):
        return mpi
    if isinstance(mpi, str):
        if mpi == "auto":
            if hasattr(os, "sched_getaffinity") and len(os.sched_getaffinity(0)) == 1:
                return None
            return find_default_mpi_runner()
        return [mpi]
    if isinstance(mpi, int):
        return find_default_mpi_runner() + ["-np", str(mpi)]
    return mpi


class Configuration(RootConfigurationContainer):
    def _generic_info(self):
        return "ASE2SPRKKR Configuration"

    def store_value_permanent(self, name: str, value, doc=None, doc_regex: bool or str = False):
        """
        Set/remove permanently the value with given name in the config file.
        No name/value checking. If the value is present, it is
        changed, otherwise, new line with setting the value is added to the file.
        If the value is to be removed, all the lines value = .... will be removed
        from the config file.

        Parameters
        ----------
        name
          Name of the value

        value
          Value to be set

        doc
          If it is not None:
           - if the value is to be removed and it is preceeded by the line containing the given
             docstring, prefixed by ' # ', the preceding line with docstring is removed as well.
           - if the value is to be added, the docstring is added before the value as well.

        doc_regex
          If True, the doc parameter can be given as regular expression.
        """
        global user_preference_file
        file = user_preferences_file()
        if not os.path.isfile(file):
            raise ValueError("Please, generate the user prefernce file using 'ase2sprkkr config -d' first.")

        with open(file, "r+") as f:
            content = f.read()
            pattern = f"(#?\\s*)*config.{name}\\s+="
            if value is None:
                line = ""
                cnt = 0
                last = ""
                if doc:
                    if not doc_regex:
                        pre = re.escape(doc)
                    else:
                        pre = doc
                    pre = f"(?:(?:^|\n)# {pre} *)?"
                else:
                    pre = ""
            else:
                line = f"config.{name} = {value}\n"
                cnt = 1
                pre = ""
                last = f"(?!(.*\n)*{pattern})"
            content, replaced = re.subn(f"{pre}(^|\n){pattern}[^\n]*(\n|$){last}", r"\1" + line, content, cnt)
            if replaced:
                f.seek(0)
                f.write(content)
            elif value is None:
                return
            else:
                i = len(content) - 1
                while i >= 0:
                    if content[i] not in ("\r", "\n", " ", "\t"):
                        break
                    i -= 1
                f.seek(i + 1)
                if i > 0:
                    f.write("\n\n")
                if doc:
                    f.write(f"# {doc}\n")
                f.write(line)
            f.truncate()


class ConfigurationOption(Option):
    def set_permanent(self, value, doc=None, doc_regex: Union[bool, str] = False):
        """Set the value and store it in the config file

        See :meth:`Configuration.store_value_permament` for the meaning of the parameters
        """

        self.set(value)
        self.store_permanent()

    def store_permanent(self, doc=None, doc_regex: Union[bool, str] = False):
        """Store the actual value in the config file

        See :meth:`Configuration.store_value_permament` for the meaning of the parameters
        """
        val = self()
        if val is not None:
            val = self._definition.type.string(val)
        config.store_value_permanent(self.get_path(), val)


class ConfigValueDefinition(ValueDefinition):
    result_class = ConfigurationOption


class GeneratedConfigValueDefinition(GeneratedValueDefinition):
    result_class = ConfigurationOption


class ConfigFileDefinition(ConfigurationRootDefinition):
    dir_common_attributes = False
    result_class = Configuration


V = ConfigValueDefinition
G = GeneratedConfigValueDefinition
MPI = CustomMixed(Bool, Array(AlwaysQString.I), Integer.I)

""" The definition of ASE2SPRKKR configuration """
definition = ConfigFileDefinition(
    "config",
    [
        V(
            "default_profile",
            AlwaysQString.I,
            is_optional=True,
            info="Explicit default execution profile. If unset, ASE and legacy profiles may be used as fallbacks.",
        ),
        ProfilesSection(
            "profiles",
            [
                V(
                    "mpi",
                    MPI,
                    is_optional=True,
                    info='MPI runner for this profile. It accepts the same values as the legacy running.mpi option.',
                ),
                V(
                    "executable_dir",
                    AlwaysQString.I,
                    is_optional=True,
                    info="Directory containing SPR-KKR executables for this profile.",
                ),
                V(
                    "executable_suffix",
                    AlwaysQString.I,
                    is_optional=True,
                    info="Suffix appended to SPR-KKR executable names for this profile.",
                ),
            ],
            is_repeated="DICT_SECTION",
            is_optional=True,
            info="Named SPR-KKR execution profiles. The reserved 'ase' and 'legacy' profiles are virtual and read-only.",
        ),
        Section(
            "running",
            [
                V(
                    "empty_spheres",
                    Keyword(
                        {
                            True: "Always do empty spheres finding.",
                            False: "Newer do empty spheres finding.",
                            "auto": "Do empty spheres finding for unconverged potential not containing any vaccuum atom.",
                        },
                        transform=None,
                        quote='"',
                    ),
                    default_value="auto",
                    info="Run empty spheres finding before calculation? Default value ``auto`` means only for SCF calculations not containing any vacuum atom.",
                ),
                V(
                    "print_output",
                    Keyword(
                        {
                            True: "Print all output of SPRKKR executables to screen.",
                            False: "Do not print any output of SPRKKR executables.",
                            "info": "Print only brief information about iterations of SCF cycle.",
                        },
                        transform=None,
                        quote='"',
                    ),
                    default_value="info",
                    info="Print output of SPRKKR calculation. Default value ``info`` prints only short info each iteration.",
                ),
                G(
                    "mpi",
                    _legacy_mpi,
                    _set_legacy_mpi,
                    is_optional=True,
                    info='Deprecated compatibility option backed by profiles["legacy"].',
                ),
                V("mpi_warning", True, info="Warn, if no MPI is found."),
            ],
            info="Default values for SPRKKR calculator parameters.",
        ),
        Section(
            "executables",
            [
                G(
                    "suffix",
                    _legacy_executable_suffix,
                    _set_legacy_executable_suffix,
                    is_optional=True,
                    info='Deprecated compatibility option backed by profiles["legacy"].',
                ),
                G(
                    "dir",
                    _legacy_executable_dir,
                    _set_legacy_executable_dir,
                    is_optional=True,
                    info='Deprecated compatibility option backed by profiles["legacy"].',
                ),
            ],
            info="Deprecated compatibility view of the legacy execution profile.",
        ),
        Section("nomad", [V("token", AlwaysQString.I, info="Token for NOMAD upload", is_optional=True)]),
    ],
)

config = definition.create_object()
