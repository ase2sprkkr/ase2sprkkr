"""Preparation and execution boundaries for SPR-KKR calculations.

``SPRKKR`` is the ASE-facing facade and owns user defaults.  This module owns
materialization of one concrete calculation and the SPR-KKR-specific adaptation
from that prepared calculation to the generic execution subsystem.
"""

from __future__ import annotations

import collections.abc
import os
import re
import tempfile
from contextlib import ExitStack
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from ..ase.symbols import filename_from_symbols
from ..bindings.empty_spheres import add_empty_spheres
from ..common.file_utils import Directory
from ..common.temporary import temporary_option_value
from ..execution.base import ExecutionResources, ResourceOwnedRunner
from ..execution.persistence import resolve_project_manager
from ..execution.postprocessing import PostprocessingRunner
from ..input_parameters.input_parameters import InputParameters
from ..potentials.potentials import Potential
from .profile import SPRKKRProfile, executor_from_profile
from .sprkkr_atoms import SPRKKRAtoms


class FilenameTemplator:
    """Render filename templates consistently for one prepared calculation.

    Placeholder values are evaluated lazily and cached, so all filenames
    created by one instance use the same date, task, symbols and counter value.
    """

    def __init__(self, label, counter, task, symbols):
        self.data = {}
        self.label = label
        self.counter = counter
        self.task = task
        self.symbols = symbols

    def _get(self, name, default):
        if name not in self.data:
            self.data[name] = default(self)
        return self.data[name]

    def _advance_counter(self):
        self.counter += 1
        return str(self.counter)

    replacements = {
        "%d": lambda self: datetime.today().strftime("%Y-%m-%d_%H:%M"),
        "%t": lambda self: str(self.task()),
        "%T": lambda self: str(self.task()).upper(),
        "%a": lambda self: str(self.symbols()),
        "%l": lambda self: self.label or "",
        "%c": lambda self: self._advance_counter(),
    }

    def __call__(self, template):
        for placeholder, value in self.replacements.items():
            if placeholder in template:
                template = template.replace(
                    placeholder, self._get(placeholder, value)
                )
        template = re.sub("^%_", "", template)
        template = re.sub("%_(%_)*", "_", template)
        return template.replace("%_", "_")


class _TemporaryNamedFile:
    """Path-backed temporary file whose pathname lives until ``close()``.

    ``NamedTemporaryFile(delete=True)`` cannot be reopened reliably by another
    process on all supported platforms. SPR-KKR and remote execution both need
    to reopen prepared files by pathname, therefore the underlying file is
    created with ``delete=False`` and this small owner removes it explicitly.
    ``close_file()`` closes only the descriptor while retaining the pathname;
    it is used after POTFIL serialization while the prepared calculation still
    owns the temporary artifact.
    """

    def __init__(self, *, directory, mode):
        self._file = tempfile.NamedTemporaryFile(
            mode=mode,
            delete=False,
            dir=directory,
        )
        self.name = self._file.name
        self._removed = False

    def __getattr__(self, name):
        return getattr(self._file, name)

    @property
    def closed(self):
        return self._file.closed

    def close_file(self):
        if not self._file.closed:
            self._file.close()

    def close(self):
        self.close_file()
        if self._removed:
            return
        self._removed = True
        try:
            os.unlink(self.name)
        except FileNotFoundError:
            pass

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass


def _calculation_path(path, directory):
    """Resolve a user path relative to the concrete calculation directory."""
    path = os.fspath(path)
    if not os.path.isabs(path):
        path = os.path.join(directory, path)
    return os.path.abspath(path)


def _named_file_name(value):
    """Return a path-backed file name, or ``None`` for anonymous streams."""
    if value is None:
        return None
    name = getattr(value, "name", value if isinstance(value, os.PathLike) else None)
    if name is None or isinstance(name, int):
        return None
    try:
        return os.fspath(name)
    except TypeError:
        return None


def _require_named_file(value, role):
    name = _named_file_name(value)
    if name is None:
        raise ValueError(
            f"SPR-KKR {role} must be backed by a filesystem path. "
            f"Pass a pathname or {role}=False to request a managed temporary file."
        )
    return name


def _open_file(
    filename,
    directory,
    templator=None,
    mode="w+",
    allow_temporary=True,
    create_subdirs=False,
):
    """Open one path-backed calculation file according to the public contract."""
    if filename is False:
        if not allow_temporary:
            raise ValueError(
                "Creation of temporary files is not allowed in this context"
            )
        return _TemporaryNamedFile(directory=directory, mode=mode)

    if not isinstance(filename, (str, os.PathLike)):
        return filename

    filename = os.fspath(filename)
    if templator:
        filename = templator(filename)
    if not os.path.isabs(filename) and directory:
        filename = os.path.join(directory, filename)
    filename = os.path.abspath(filename)
    if create_subdirs:
        parent = os.path.dirname(filename)
        if parent:
            os.makedirs(parent, exist_ok=True)
    return open(filename, mode)


def _own_file(stack, owned_ids, file):
    """Register ownership of a file-like object exactly once."""
    if file is None or id(file) in owned_ids:
        return
    close = getattr(file, "close", None)
    if close is None:
        return
    owned_ids.add(id(file))
    stack.callback(close)


def _resolve_potential_and_atoms(
    *,
    potential,
    atoms,
    potential_file,
    default_potential,
    default_atoms,
    restart_scf,
    directory,
):
    """Resolve the concrete structure/potential pair for preparation."""
    save_input = None
    if potential is None and atoms is None:
        potential = default_potential
    if potential is None:
        potential = atoms is not None or default_atoms is not None

    potential_specified = potential is not None and potential is not True
    if atoms is not None and potential_specified:
        raise ValueError(
            "You can not provide both a potential and atoms object to the SPRKKR "
            "calculate method; one is generated from the other."
        )

    atoms = default_atoms if atoms is None else atoms
    if potential is True:
        if atoms is None:
            raise ValueError(
                "Potential is True, which requests generation from atoms, but no "
                "atoms object has been supplied."
            )
        potential = Potential.from_atoms(atoms)
    elif potential is False:
        if potential_file is not None:
            raise ValueError(
                "When potential=False, POTFIL from the input parameters is used; "
                "potential_file therefore must not be specified."
            )
        potential = None
    elif potential is None:
        raise ValueError(
            "No potential is available. Supply a Potential, atoms, potential=True "
            "to generate one from atoms, or potential=False to use POTFIL from the "
            "input parameters."
        )
    elif isinstance(potential, (str, os.PathLike)):
        source = _calculation_path(potential, directory)
        if potential_file is not None or restart_scf:
            potential = Potential.from_file(source)
        else:
            potential_file = source
            potential = False
            save_input = "check_potential"

    if isinstance(potential, Potential):
        atoms = potential.atoms
    elif atoms is not None:
        atoms = SPRKKRAtoms.promote_ase_atoms(atoms)

    return potential, atoms, potential_file, save_input


def _potential_from_input_file(parameters, parameter_filename):
    potfil = parameters.CONTROL.POTFIL.result
    if not potfil:
        return None
    return _calculation_path(potfil, os.path.dirname(parameter_filename))


def _resolve_input_parameters(
    *,
    input_parameters,
    default_input_parameters,
    input_file,
    potential_file,
    potential,
    save_input,
    options,
    task,
    directory,
):
    """Resolve and snapshot task parameters for one concrete calculation."""
    value = default_input_parameters if input_parameters is None else input_parameters
    potential_missing = potential is None or potential is False

    if isinstance(value, (str, os.PathLike)):
        if isinstance(value, str) and InputParameters.is_task_name(value):
            value = InputParameters.create(value)
        else:
            parameter_filename = _calculation_path(value, directory)
            parsed = InputParameters.from_file(parameter_filename)
            source_potential = _potential_from_input_file(parsed, parameter_filename)

            if save_input == "check_potential":
                save_input = source_potential != potential_file
            elif (
                save_input is None
                and potential_missing
                and potential_file is None
                and not options
                and not task
                and input_file is None
            ):
                save_input = False

            if not save_input and (options or task or input_file is not None):
                save_input = True

            if save_input is False:
                input_file = parameter_filename
            elif potential_missing and potential_file is None and source_potential:
                potential_file = source_potential
            value = parsed

    if value is None:
        value = default_input_parameters
    if value is None:
        raise ValueError("No input parameters are available for the calculation")

    input_parameters = value.copy(copy_values=True)
    if task:
        input_parameters.change_task(task)
    if options:
        input_parameters.set(options, unknown="find")

    if (
        not input_parameters.CONTROL.POTFIL.result
        and potential_missing
        and not potential_file
    ):
        raise ValueError(
            "Potential in the input parameters is not set and no atoms or "
            "Potential object has been supplied."
        )

    if save_input is None or save_input == "check_potential":
        save_input = True
    return input_parameters, input_file, potential_file, bool(save_input)


def _resolve_empty_spheres(
    *, input_parameters, atoms, empty_spheres, default_empty_spheres
):
    """Apply the configured empty-sphere preparation policy."""
    if empty_spheres is None:
        empty_spheres = default_empty_spheres
    if empty_spheres == "auto":
        empty_spheres = (
            input_parameters.task_name == "SCF"
            and atoms is not None
            and (
                not atoms.has_potential()
                or atoms.potential.SCF_INFO.SCFSTATUS() == "START"
            )
        )
        if empty_spheres:
            empty_spheres = not any(site.is_vacuum() for site in atoms.sites)

    if not empty_spheres or atoms is None:
        return

    if isinstance(empty_spheres, tuple):
        if len(empty_spheres) != 2:
            raise ValueError(
                "empty_spheres tuple must contain (min_radius, max_radius)"
            )
        empty_spheres = {
            "min_radius": empty_spheres[0],
            "max_radius": empty_spheres[1],
        }
    if not isinstance(empty_spheres, collections.abc.Mapping):
        empty_spheres = {}

    add_empty_spheres(atoms, **empty_spheres)


def _symbols_for_filenames(
    atoms,
    input_parameters,
    input_file,
    potential_file,
    default_input_file,
    default_potential_file,
):
    """Return a stable best-effort system name for filename templates."""
    if atoms is not None:
        return filename_from_symbols(atoms.symbols)
    for name in "POTFIL", "DATASET":
        option = input_parameters.CONTROL[name]
        if option.is_set():
            return Path(option()).stem
    for candidate in (
        input_file,
        potential_file,
        default_input_file,
        default_potential_file,
    ):
        if isinstance(candidate, (str, os.PathLike)):
            stem = Path(candidate).stem
            if "%" not in stem:
                return stem
    return "sprkkr"


@dataclass
class PreparedCalculation:
    """One fully materialized SPR-KKR calculation.

    Instances own the concrete parameter snapshot, path-backed calculation
    files and, when applicable, an entered temporary working-directory context.
    Ownership can be transferred to a runner, explicitly released for the
    legacy ``save_input(return_files=True)`` API, or terminated with
    :meth:`close`.

    :meth:`prepare` is the canonical factory because it resolves all public
    calculator sentinels and establishes these invariants. The generated
    dataclass constructor validates path-backed file invariants as well, so it
    cannot silently create an execution-incompatible instance.
    """

    input_parameters: InputParameters
    input_file: object
    potential_file: str | None
    output_file: object | None
    directory: str
    counter: int
    _file_resources: ExitStack | None
    _directory_resources: ExitStack | None
    _temporary_potential: bool = False
    _closed: bool = False

    def __post_init__(self):
        _require_named_file(self.input_file, "input_file")
        if self.output_file is not None:
            _require_named_file(self.output_file, "output_file")
        if self.potential_file is not None:
            try:
                os.fspath(self.potential_file)
            except TypeError:
                raise ValueError(
                    "SPR-KKR potential_file must be backed by a filesystem path"
                ) from None
        self.directory = os.path.abspath(os.fspath(self.directory))

    @classmethod
    def prepare(
        cls,
        *,
        default_atoms,
        default_potential,
        default_input_parameters,
        default_input_file,
        default_potential_file,
        default_output_file,
        default_directory,
        default_empty_spheres,
        label,
        counter,
        atoms=None,
        input_parameters=None,
        potential=None,
        input_file=None,
        potential_file=None,
        output_file=False,
        directory=None,
        create_subdirs=False,
        options=None,
        task=None,
        empty_spheres=None,
        allow_temporary_files=False,
        restart_scf=False,
    ) -> PreparedCalculation:
        """Materialize all state and files required for one SPR-KKR run.

        ``default_*`` arguments are supplied by the calculator facade; explicit
        per-call arguments override them. On success the returned object owns
        every resource opened by preparation. On failure every resource
        acquired up to that point is released before the original exception is
        propagated.
        """
        options = options or {}
        directory_resources = ExitStack()
        file_resources = ExitStack()
        owned_file_ids = set()
        temporary_potential = False

        try:
            directory_context = Directory.new(directory, default=default_directory)
            work_directory = directory_resources.enter_context(directory_context).path
            work_directory = os.path.abspath(os.fspath(work_directory))

            potential, atoms, potential_file, save_input = (
                _resolve_potential_and_atoms(
                    potential=potential,
                    atoms=atoms,
                    potential_file=potential_file,
                    default_potential=default_potential,
                    default_atoms=default_atoms,
                    restart_scf=restart_scf,
                    directory=work_directory,
                )
            )
            (
                input_parameters,
                input_file,
                potential_file,
                save_input,
            ) = _resolve_input_parameters(
                input_parameters=input_parameters,
                default_input_parameters=default_input_parameters,
                input_file=input_file,
                potential_file=potential_file,
                potential=potential,
                save_input=save_input,
                options=options,
                task=task,
                directory=work_directory,
            )

            task_definition = input_parameters.task_definition_required()
            task_name = task_definition.name
            if restart_scf:
                if task_name != "SCF":
                    raise ValueError(
                        "restart_scf=True is only valid for an SCF calculation."
                    )
                if not isinstance(potential, Potential):
                    raise ValueError(
                        "restart_scf=True requires a writable Potential or atoms "
                        "from which one can be generated."
                    )

            _resolve_empty_spheres(
                input_parameters=input_parameters,
                atoms=atoms,
                empty_spheres=empty_spheres,
                default_empty_spheres=default_empty_spheres,
            )

            scf_status_override = None
            if (
                isinstance(potential, Potential)
                and task_name != "SCF"
                and potential.SCF_INFO.SCFSTATUS() == "START"
                and potential.is_scf_converged(atoms)
            ):
                scf_status_override = "CONVERGED"

            templator = FilenameTemplator(
                label,
                counter,
                lambda: task_name,
                lambda: _symbols_for_filenames(
                    atoms,
                    input_parameters,
                    input_file,
                    potential_file,
                    default_input_file,
                    default_potential_file,
                ),
            )

            selected_input = default_input_file if input_file is None else input_file
            if selected_input is True:
                selected_input = "%l%_%a_%T.inp"
            opened_input = _open_file(
                selected_input,
                work_directory,
                templator,
                allow_temporary=allow_temporary_files,
                create_subdirs=create_subdirs,
                mode="w+" if save_input else "r",
            )
            _require_named_file(opened_input, "input_file")
            _own_file(file_resources, owned_file_ids, opened_input)

            def from_input_name(template, replacement, default):
                if template is not True:
                    return template
                if not save_input:
                    name = _named_file_name(opened_input)
                    if name:
                        base = os.path.basename(name)
                        if base.endswith(".inp"):
                            return base[:-4] + replacement
                        if base.endswith(".in"):
                            return base[:-3] + replacement
                        return base + replacement
                return default

            if potential is not None and potential is not False:
                selected_potential = (
                    default_potential_file
                    if potential_file is None
                    else potential_file
                )
                selected_potential = from_input_name(
                    selected_potential, ".pot", "%a.pot"
                )
                temporary_potential = selected_potential is False
                potential_handle = _open_file(
                    selected_potential,
                    work_directory,
                    templator,
                    allow_temporary=allow_temporary_files,
                    create_subdirs=create_subdirs,
                )
                potential_name = _require_named_file(
                    potential_handle, "potential_file"
                )
                _own_file(file_resources, owned_file_ids, potential_handle)

                if scf_status_override is None:
                    potential.save_to_file(
                        potential_handle, atoms, restart_scf=restart_scf
                    )
                else:
                    with temporary_option_value(
                        potential.SCF_INFO.SCFSTATUS, scf_status_override
                    ):
                        potential.save_to_file(
                            potential_handle, atoms, restart_scf=restart_scf
                        )
                flush = getattr(potential_handle, "flush", None)
                if flush is not None:
                    flush()
                if isinstance(potential_handle, _TemporaryNamedFile):
                    potential_handle.close_file()
                else:
                    potential_handle.close()
                potential_file = potential_name

            if save_input:
                if potential_file:
                    potential_file = os.path.abspath(os.fspath(potential_file))
                    input_parameters.CONTROL.POTFIL = potential_file
                    input_parameters.CONTROL.POTFIL.result = os.path.relpath(
                        potential_file, work_directory
                    )
                if not input_parameters.CONTROL.DATASET.is_set():
                    input_name = _named_file_name(opened_input)
                    if input_name:
                        input_parameters.CONTROL.DATASET.result = Path(
                            input_name
                        ).stem
                input_parameters.save_to_file(opened_input, atoms)
                opened_input.seek(0)
            else:
                input_name = _require_named_file(opened_input, "input_file")
                input_dir = os.path.dirname(input_name) or work_directory
                potential_file = _calculation_path(
                    input_parameters.CONTROL.POTFIL.result,
                    input_dir,
                )

            opened_output = None
            selected_output = (
                default_output_file if output_file is None else output_file
            )
            if selected_output is not False and selected_output is not None:
                selected_output = from_input_name(
                    selected_output, ".out", "%a_%T.out"
                )
                opened_output = _open_file(
                    selected_output,
                    work_directory,
                    templator,
                    mode="wb",
                    allow_temporary=allow_temporary_files,
                    create_subdirs=create_subdirs,
                )
                _require_named_file(opened_output, "output_file")
                _own_file(file_resources, owned_file_ids, opened_output)

            return cls(
                input_parameters=input_parameters,
                input_file=opened_input,
                potential_file=potential_file,
                output_file=opened_output,
                directory=work_directory,
                counter=templator.counter,
                _file_resources=file_resources,
                _directory_resources=directory_resources,
                _temporary_potential=temporary_potential,
            )
        except BaseException:
            try:
                file_resources.close()
            except BaseException:
                pass
            try:
                directory_resources.close()
            except BaseException:
                pass
            raise

    def create_run(
        self,
        *,
        result_context=None,
        profile: SPRKKRProfile,
        executor=None,
        project=None,
        print_output=False,
        read_callback=None,
        callback=None,
        gdb=False,
        error_file=None,
        create_subdirs=False,
    ):
        """Create an unstarted runner and transfer this calculation's resources.

        The prepared calculation does not depend on the calculator facade.
        ``result_context`` is optional runtime context for output readers; the
        :class:`SPRKKR` facade supplies itself so result convenience APIs can
        construct follow-up calculations, while standalone prepared runs may
        leave the context unset.

        ``error_file=None`` stores stderr as ``<output_file>.err``. ``True``
        explicitly merges stderr into ``output_file``; ``False`` requests a
        managed temporary stderr file; a path or file object is used directly.
        """
        return _create_sprkkr_run(
            self,
            result_context=result_context,
            profile=profile,
            executor=executor,
            project=project,
            print_output=print_output,
            read_callback=read_callback,
            callback=callback,
            gdb=gdb,
            error_file=error_file,
            create_subdirs=create_subdirs,
        )

    @property
    def closed(self):
        """Whether this object has relinquished or closed all owned resources."""
        return self._closed

    @property
    def temporary_potential(self):
        """Whether POTFIL is a managed temporary pathname."""
        return self._temporary_potential

    def release_file_ownership(self):
        """Transfer ownership of open input/output files to the caller.

        The files themselves are not closed. This supports the legacy
        ``save_input(return_files=True)`` API. A temporary POTFIL cannot be
        transferred by that API because it returns only its pathname; callers
        must supply a persistent POTFIL name in that case.
        """
        stack = self._file_resources
        self._file_resources = None
        if stack is not None:
            stack.pop_all()
        self._update_closed_state()

    def close_directory_context(self):
        """Close the entered working-directory context after file transfer."""
        stack = self._directory_resources
        self._directory_resources = None
        if stack is not None:
            stack.close()
        self._update_closed_state()

    def _update_closed_state(self):
        if self._file_resources is None and self._directory_resources is None:
            self._closed = True

    def close(self):
        """Close all resources still owned by this prepared calculation."""
        if self._closed:
            return

        file_resources = self._file_resources
        directory_resources = self._directory_resources
        self._file_resources = None
        self._directory_resources = None
        self._closed = True

        error = None
        if file_resources is not None:
            try:
                file_resources.close()
            except BaseException as exc:
                error = exc
        if directory_resources is not None:
            try:
                directory_resources.close()
            except BaseException as exc:
                if error is None:
                    error = exc
        if error is not None:
            raise error


def _prepare_error_file(error_file, output_file, directory, create_subdirs):
    """Materialize stderr for one run without changing preparation semantics."""
    if error_file is True or error_file is output_file:
        return output_file, False

    output_name = _named_file_name(output_file)
    if error_file is None:
        if output_name is None:
            return None, False
        error_file = output_name + ".err"
    elif isinstance(error_file, (str, os.PathLike)) and output_name is not None:
        error_name = _calculation_path(error_file, directory)
        if os.path.abspath(output_name) == error_name:
            return output_file, False

    opened_error = _open_file(
        error_file,
        directory,
        mode="wb",
        allow_temporary=True,
        create_subdirs=create_subdirs,
    )
    _require_named_file(opened_error, "error_file")
    return opened_error, opened_error is not output_file


def _create_sprkkr_run(
    prepared: PreparedCalculation,
    *,
    result_context,
    profile: SPRKKRProfile,
    executor=None,
    project=None,
    print_output=False,
    read_callback=None,
    callback=None,
    gdb=False,
    error_file=None,
    create_subdirs=False,
):
    """Adapt a prepared SPR-KKR calculation to the execution subsystem.

    This SPR-KKR layer owns task metadata, reader construction and TaskResult
    reconstruction. Executors see only :class:`ExecutionSpec` and produce raw
    :class:`ExecutionOutcome` values; the generic postprocessing layer merely
    invokes the callable defined here.
    """
    try:
        opened_error, owns_error = _prepare_error_file(
            error_file,
            prepared.output_file,
            prepared.directory,
            create_subdirs,
        )
    except BaseException:
        prepared.close()
        raise

    resources = ExecutionResources(
        prepared,
        opened_error if owns_error else None,
    )
    handed_off = False
    input_parameters = prepared.input_parameters

    try:
        task = input_parameters.task_definition_required()
        reader = task.output_reader(
            input_parameters=input_parameters,
            calculator=result_context,
            print_output=print_output,
            directory=prepared.directory,
        )
        if executor is None:
            executor = executor_from_profile(profile)
        manager = resolve_project_manager(project)

        execution_profile = profile
        if task.mpi:
            execution_profile = profile.with_overrides(
                mpi=executor.resolve_mpi(profile.mpi)
            )

        input_files = tuple(
            file
            for file in (
                prepared.input_file,
                prepared.potential_file,
                *input_parameters.input_files(prepared.directory),
            )
            if file is not None
        )
        output_files = tuple(
            file
            for file in (
                prepared.output_file,
                opened_error if opened_error is not prepared.output_file else None,
                *input_parameters.output_files(prepared.directory),
            )
            if file is not None
        )
        cache_files = tuple(
            file
            for file in input_parameters.cache_files(prepared.directory)
            if file is not None
        )

        spec = execution_profile.create_execution(
            task.executable,
            prepared.input_file,
            prepared.output_file,
            prepared.directory,
            error_file=opened_error,
            input_files=input_files,
            output_files=output_files,
            cache_files=cache_files,
            mpi_supported=task.mpi,
            gdb=gdb,
        )

        def process_result(outcome):
            return reader.read_from_file(
                outcome.output_file,
                outcome.error_file,
                return_code=(
                    0 if outcome.return_code is None else outcome.return_code
                ),
                input_file=outcome.input_file,
            )

        def finalized(out):
            resources.close(suppress_errors=True)
            if callback:
                callback(out)

        if manager is None:
            raw_runner = executor.create_run(spec)
            runner = PostprocessingRunner(
                raw_runner,
                process_result,
                callback=finalized,
            )
        else:
            runner = manager.create_run(
                executor,
                spec,
                process_result=process_result,
                task=task.name,
                input_definition=input_parameters.definition_name,
                dataset=input_parameters.CONTROL.DATASET.result,
                callback=finalized,
            )
        if read_callback:
            runner.add_output_listener(read_callback)

        runner = ResourceOwnedRunner(runner, resources)
        handed_off = True
        return runner
    finally:
        if not handed_off:
            resources.close(suppress_errors=True)
