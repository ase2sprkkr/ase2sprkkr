"""ASE2SPRKKR user configuration."""

import os
from collections.abc import Mapping
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
    section._container.profiles._set_legacy_value(name, value)


def _legacy_mpi(section):
    return _legacy_value(section, "mpi", "auto")


def _set_legacy_mpi(section, value):
    _set_legacy_value(section, value, "mpi")


def _legacy_mpi_warning(section):
    return _legacy_value(section, "mpi_warning", True)


def _set_legacy_mpi_warning(section, value):
    _set_legacy_value(section, value, "mpi_warning")


def _legacy_executable_dir(section):
    return _legacy_value(section, "executable_dir")


def _set_legacy_executable_dir(section, value):
    _set_legacy_value(section, value, "executable_dir")


def _legacy_executable_suffix(section):
    return _legacy_value(section, "executable_suffix", _get_suffix)


def _set_legacy_executable_suffix(section, value):
    _set_legacy_value(section, value, "executable_suffix")


def user_preferences_directory():
    """Return the ASE2SPRKKR user configuration directory."""
    return platformdirs.user_config_dir("ase2sprkkr", "ase2sprkkr")


def user_preferences_file():
    """Return the preferred YAML user configuration filename."""
    return os.path.join(user_preferences_directory(), "config.yaml")


def legacy_user_preferences_file():
    """Return the optional Python user-configuration hook filename."""
    return os.path.join(user_preferences_directory(), "__init__.py")


def _load_legacy_user_preferences(file):
    import types
    import importlib.machinery

    loader = importlib.machinery.SourceFileLoader("ase2sprkkr.personal", file)
    mod = types.ModuleType(loader.name)
    mod.config = config
    loader.exec_module(mod)


def _profile_overlay_values(yaml_file):
    """Return profile values from the Python hook not overridden by YAML."""
    from .configuration_yaml import load_yaml_document

    python_profiles = {
        name: profile.as_dict(only_changed="explicit", generated=False) or {}
        for name, profile in config.profiles._raw_items()
    }
    if not python_profiles:
        return {}

    document = load_yaml_document(yaml_file)
    yaml_profiles = document.get("profiles", {})
    if not isinstance(yaml_profiles, Mapping):
        return {}

    overlay = {}
    for name, values in python_profiles.items():
        if name not in yaml_profiles:
            overlay[name] = values
            continue
        yaml_values = yaml_profiles[name]
        if not isinstance(yaml_values, Mapping):
            continue
        missing = {key: value for key, value in values.items() if key not in yaml_values}
        if missing:
            overlay[name] = missing
    return overlay


def load_user_preferences():
    """Load the Python hook and then overlay ``config.yaml`` if present.

    Expected invalid YAML values are handled by ``load_configuration_yaml`` at
    the narrow boundary where they can be recovered safely. Unexpected loader
    and Python-hook exceptions deliberately propagate instead of silently
    falling back to a partially/default configured process.
    """
    yaml_file = user_preferences_file()
    legacy_file = legacy_user_preferences_file()

    if os.path.isfile(legacy_file):
        _load_legacy_user_preferences(legacy_file)
    if os.path.isfile(yaml_file):
        from .configuration_yaml import load_configuration_yaml

        overlay = _profile_overlay_values(yaml_file)
        load_configuration_yaml(config, yaml_file)
        if overlay:
            with config.profiles._mutation("set") as (transaction, _policy):
                config.profiles._stage_persisted(
                    transaction,
                    overlay,
                    unknown="fail",
                    merge=True,
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
    if config.running.mpi_warning():
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

    def save_user_preferences(self):
        """Persist explicit configuration values to the YAML user config."""
        from .configuration_yaml import load_yaml_document, save_configuration_yaml

        file = user_preferences_file()
        document = load_yaml_document(file) if os.path.isfile(file) else None
        save_configuration_yaml(self, file, document)


class ConfigurationOption(Option):
    def set_permanent(self, value, doc=None, doc_regex=False):
        """Set the value and persist the complete user configuration."""
        self.set(value)
        self.store_permanent()

    def store_permanent(self, doc=None, doc_regex=False):
        """Persist the current configuration to the YAML user config."""
        self._get_root_container().save_user_preferences()


class GeneratedConfigurationOption(ConfigurationOption):
    """Configuration option whose value is backed by another config object."""


class ConfigValueDefinition(ValueDefinition):
    result_class = ConfigurationOption


class GeneratedConfigValueDefinition(GeneratedValueDefinition):
    result_class = GeneratedConfigurationOption


class ConfigFileDefinition(ConfigurationRootDefinition):
    dir_common_attributes = False
    result_class = Configuration


V = ConfigValueDefinition
G = GeneratedConfigValueDefinition
MPI = CustomMixed(
    (Bool.I, Array(AlwaysQString.I), Integer.I, AlwaysQString.I),
    string_type=AlwaysQString.I,
)
SYNC_POLICY = Keyword(
    {
        "changed": "Transfer files only when their content differs.",
        "always": "Always transfer files.",
        "missing": "Transfer files only when the destination is missing.",
        "never": "Never transfer files.",
    },
    transform=None,
)
EXECUTOR = Keyword(
    {
        "local": "Run directly on this machine.",
        "ssh": "Run through ASE2SPRKKR on a machine reached by SSH.",
        "psij": "Submit through a PSI/J scheduler backend.",
    },
    transform=None,
)

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
                    default_value="auto",
                    is_optional=True,
                    info='MPI runner for this profile. Defaults to "auto" and accepts the same values as the legacy running.mpi option.',
                ),
                V(
                    "mpi_warning",
                    True,
                    info="Whether missing MPI support should be reported for this profile.",
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
                V(
                    "executor",
                    EXECUTOR,
                    is_optional=True,
                    info="Execution backend for this profile. If unset, local execution is used.",
                ),
                V(
                    "host",
                    AlwaysQString.I,
                    is_optional=True,
                    info="SSH destination for an ssh executor, for example user@cluster.",
                ),
                V(
                    "local_root",
                    AlwaysQString.I,
                    is_optional=True,
                    info="Local root mapped into the remote workspace by an ssh executor.",
                ),
                V(
                    "remote_root",
                    AlwaysQString.I,
                    is_optional=True,
                    info="Remote root used for persistent workspaces by an ssh executor.",
                ),
                V(
                    "ssh_options",
                    Array(AlwaysQString.I),
                    is_optional=True,
                    info="Additional arguments passed to ssh, such as ['-p', '2222'].",
                ),
                V(
                    "remote_command",
                    Array(AlwaysQString.I),
                    is_optional=True,
                    info="Command prefix used to invoke the ASE2SPRKKR remote helper.",
                ),
                V(
                    "input_sync",
                    SYNC_POLICY,
                    is_optional=True,
                    info="Input-file synchronization policy for an ssh executor.",
                ),
                V(
                    "output_sync",
                    SYNC_POLICY,
                    is_optional=True,
                    info="Output-file synchronization policy for an ssh executor.",
                ),
                V(
                    "psij_backend",
                    AlwaysQString.I,
                    is_optional=True,
                    info="PSI/J executor backend, for example slurm or pbs.",
                ),
                V(
                    "psij_url",
                    AlwaysQString.I,
                    is_optional=True,
                    info="Optional PSI/J executor URL.",
                ),
            ],
            is_repeated="DICT_SECTION",
            is_optional=True,
            info="Named SPR-KKR execution profiles. The reserved 'ase' profile is virtual; 'legacy' is persisted but read-only through the profiles API.",
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
                G(
                    "mpi_warning",
                    _legacy_mpi_warning,
                    _set_legacy_mpi_warning,
                    is_optional=True,
                    info='Deprecated compatibility option backed by profiles["legacy"].',
                ),
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
