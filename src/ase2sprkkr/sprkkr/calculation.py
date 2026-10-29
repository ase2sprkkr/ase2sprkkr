"""Prepared SPR-KKR calculations and their execution adapter.

``SPRKKR`` supplies user defaults; ``_preparation`` resolves state and writes
files. This module orchestrates preparation, owns the resulting resources,
and adapts prepared calculations to the generic execution subsystem.
"""

from __future__ import annotations

import os
from contextlib import ExitStack
from dataclasses import dataclass

from ..common.file_utils import Directory
from ..common.output_printing import output_listeners
from ..execution.base import ExecutionResources, ResourceOwnedRunner
from ..execution.persistence import resolve_project_manager
from ..execution.postprocessing import PostprocessingRunner
from ..input_parameters.input_parameters import InputParameters
from . import _preparation as preparation
from ._preparation import FilenameTemplator as FilenameTemplator
from .profile import SPRKKRProfile, executor_from_profile


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
        preparation.require_named_file(self.input_file, "input_file")
        if self.output_file is not None:
            preparation.require_named_file(self.output_file, "output_file")
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

        Relative filenames and references inside input parameters use the
        calculation directory as their base. Files opened here have absolute
        names; filename mappings serialize references relative to that
        directory, including when the input file is stored in a subdirectory.
        """
        directory_resources = ExitStack()
        file_resources = ExitStack()

        try:
            directory_context = Directory.new(directory, default=default_directory)
            work_directory = directory_resources.enter_context(directory_context).path
            work_directory = os.path.abspath(os.fspath(work_directory))

            resolved_potential = preparation.resolve_potential_and_atoms(
                potential=potential,
                atoms=atoms,
                potential_file=potential_file,
                default_potential=default_potential,
                default_atoms=default_atoms,
                restart_scf=restart_scf,
                directory=work_directory,
            )
            resolved_input = preparation.resolve_input_parameters(
                input_parameters=input_parameters,
                default_input_parameters=default_input_parameters,
                input_file=input_file,
                potential=resolved_potential,
                options=options,
                task=task,
                directory=work_directory,
            )
            parameters = resolved_input.parameters
            atoms = resolved_potential.atoms
            task_name = parameters.task_definition_required().name
            preparation.validate_restart(resolved_potential.potential, task_name, restart_scf)
            preparation.resolve_empty_spheres(
                input_parameters=parameters,
                atoms=atoms,
                empty_spheres=empty_spheres,
                default_empty_spheres=default_empty_spheres,
            )

            templator = preparation.FilenameTemplator(
                label,
                counter,
                lambda: task_name,
                lambda: preparation.symbols_for_filenames(
                    atoms,
                    parameters,
                    resolved_input.filename,
                    resolved_input.potential_file,
                    default_input_file,
                    default_potential_file,
                ),
            )
            files = preparation.CalculationFiles(
                work_directory,
                templator,
                file_resources,
                allow_temporary=allow_temporary_files,
                create_subdirs=create_subdirs,
            )
            opened_input = files.open_input(resolved_input, default_input_file)
            reused_input = None if resolved_input.write else preparation.named_file_name(opened_input)

            potential_file, temporary_potential = preparation.write_potential(
                resolved_potential,
                resolved_input.potential_file,
                files,
                default_filename=default_potential_file,
                task_name=task_name,
                restart_scf=restart_scf,
            )
            preparation.write_input(resolved_input, opened_input, potential_file, atoms)

            opened_output = files.open_output(
                output_file,
                default_filename=default_output_file,
                reused_input=reused_input,
            )

            return cls(
                input_parameters=parameters,
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

    output_name = preparation.named_file_name(output_file)
    if error_file is None:
        if output_name is None:
            return None, False
        error_file = output_name + ".err"
    elif isinstance(error_file, (str, os.PathLike)) and output_name is not None:
        error_name = preparation.calculation_path(error_file, directory)
        if os.path.abspath(output_name) == error_name:
            return output_file, False

    opened_error = preparation.open_file(
        error_file,
        directory,
        mode="wb",
        allow_temporary=True,
        create_subdirs=create_subdirs,
    )
    preparation.require_named_file(opened_error, "error_file")
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
        listeners, print_info = output_listeners(print_output)
        task = input_parameters.task_definition_required()
        reader = task.output_reader(
            input_parameters=input_parameters,
            calculator=result_context,
            print_info=print_info,
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

        def flush_listeners():
            for listener in listeners:
                listener.flush()

        def process_result(outcome):
            flush_listeners()
            return reader.read_from_file(
                outcome.output_file,
                outcome.error_file,
                return_code=(
                    0 if outcome.return_code is None else outcome.return_code
                ),
                input_file=outcome.input_file,
            )

        def finalized(out):
            flush_listeners()
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
        for listener in listeners:
            runner.add_output_listener(listener)
        if read_callback:
            runner.add_output_listener(read_callback)

        runner = ResourceOwnedRunner(runner, resources)
        handed_off = True
        return runner
    finally:
        if not handed_off:
            resources.close(suppress_errors=True)
