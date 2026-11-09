Execution, persistent runs and KKRunner
=========================================

ASE2SPRKKR prepares SPR-KKR inputs and interprets its outputs. The separate
`KKRunner library <https://ase2sprkkr.github.io/kkrunner/>`_ owns local/MPI,
SSH and scheduler execution, file transfers and durable run history.
It has no dependency on ASE2SPRKKR or on its scientific readers.

For a normal calculation, continue using ``SPRKKR.calculate()``. Use
``SPRKKR.create_run()`` when you need to separate preparation/submission from
waiting, observe progress, or retain a handle for a persistent calculation.
The returned application result is a ``TaskResult``; KKRunner's raw
``ExecutionOutcome`` describes process facts before application interpretation.

Local calculations
--------------------

.. code-block:: python

   from ase.build import bulk
   from ase2sprkkr import SPRKKR

   calculator = SPRKKR(atoms=bulk("Al", "fcc", a=4.05), directory="al-scf")
   run = calculator.create_run()
   run.submit()
   result = run.wait()

Preparation is performed when ``create_run()`` is called; submission starts the
executable. The run owns the resources transferred by preparation. Retain its
working directory until execution, transfers and result processing are complete.
These examples require installed SPR-KKR executables.

Persistent project calculations
---------------------------------

Persistence is explicit. Create a project and pass a project-bound manager;
omitting ``project`` does not silently enable global tracking.

.. code-block:: python

   from pathlib import Path
   from kkrunner.projects import ProjectRepository
   from ase2sprkkr.sprkkr.projects import project_execution_manager

   projects = ProjectRepository(Path("calculation-state/projects").absolute())
   project = projects.create("Aluminium", directory=Path("al-scf").absolute())
   manager = project_execution_manager(projects, project=project)
   run = calculator.create_run(project=manager)
   run.submit()
   run_id = run.id
   result = run.wait()

The application helper constructs KKRunner's ``ProjectExecutionManager`` using
ASE2SPRKKR's storage policy. ``RunManager`` is the distinct project-free KKRunner
API for clients that supply an explicit ``RunStore``. Use those actual class
names when extending the application.

Project definitions retain ASE2SPRKKR's existing application-state location.
New-format run records use KKRunner's state namespace. Passing an explicit
project repository retains its colocated run store. Old development record
formats are not migrated or silently ignored as current records.

Reconnect and load history
----------------------------

Reopen the same project repository and its manager after restarting the client.
Use the application adapter when reconnecting, so it supplies the correct reader:

.. code-block:: python

   from ase2sprkkr.outputs.project_results import (
       attach_project_run, load_project_result, retry_project_postprocessing,
   )

   # For an active, downloading or waiting-for-processing run:
   run = attach_project_run(manager, run_id)
   result = run.wait()

   # For a terminal historical run:
   result = load_project_result(manager, run_id)

These are alternative actions according to the record's state. Do not attach a
terminal run merely to load its result. ``TaskResult.from_run(run_id,
manager=manager)`` is another historical-loading entry point.

``attach_project_run`` also reopens a supported failed download before attaching.
For failed or interrupted result processing, call
``retry_project_postprocessing(manager, run_id)`` explicitly. It verifies the
preserved files and parses them again without resubmitting SPR-KKR.

Recovery after interruption during submission depends on the selected backend.
Managed local execution uses retained service receipts. Slurm can also discover
jobs by their run marker; discovery of completed jobs needs accounting that
stores job comments. Other PSI/J backends do not implement discovery without a
saved backend ID, even when attaching to a known job works. Uncertain acceptance
leaves the record reserved in ``submitting``; it does not trigger a duplicate
SPR-KKR calculation. See the
`backend recovery capabilities <https://ase2sprkkr.github.io/kkrunner/backends.html#backend-recovery>`_
for target-side requirements and limits.

``stop_observation()`` detaches the observer; ``cancel()`` requests cancellation
of the backend. ``UNKNOWN`` means observation is uncertain, not that execution
has finished. Required processing that is unavailable or owned by another client
raises ``ProcessingPending``; callers must not treat a raw process outcome as a
replacement for a ``TaskResult``.

Only the application's namespaced metadata selects a task definition and reader.
Generic or foreign application records in a shared store are not SPR-KKR history.
Historical loading verifies captured artifacts; overwritten files must be reported
rather than parsed as the original run's results.

Profiles, remote execution and schedulers
-------------------------------------------

A profile describes the backend connection and per-run launch/staging policy.
An explicit per-call executor takes precedence. An explicit per-call profile
selects its backend; otherwise a calculator-level executor takes precedence over
the calculator's default profile backend.

.. code-block:: python

   from ase2sprkkr.sprkkr.profile import SPRKKRProfile

   profile = SPRKKRProfile(
       executor="ssh",
       host="user@cluster",
       local_root="/absolute/local/work",
       remote_root="/absolute/remote/work",
       mpi="auto",
   )
   run = calculator.create_run(profile=profile, project=manager)

Install the same KKRunner release on the target and make ``kkrunner remote``
available to noninteractive SSH commands. The target needs SPR-KKR and its MPI
runtime when requested; it does not need ASE2SPRKKR merely to execute its command.
For remote scheduler execution, set ``psij_backend="slurm"`` and install
``kkrunner[psij]`` on the target. ``remote_command`` can explicitly select an
absolute command such as ``["/opt/kkrunner/bin/kkrunner", "remote"]``. It must
not select the end-user
``ase2sprkkr remote`` subset as the execution protocol.

Use a consistent ``filesystem_id`` in profiles that address the same filesystem
through different host aliases or shared mounts. KKRunner cannot infer that
relationship from hostnames. Required MPI cannot fall back to serial execution;
auto-detection may allow fallback according to the selected MPI policy.

Files, manifests and resource lifetime
----------------------------------------

ASE2SPRKKR resolves task declarations into KKRunner's ``TaskFile`` and manifest
data. Scientific filename templates, parser registries and format descriptions
belong to ASE2SPRKKR. ``TaskFile.type`` and ``info`` are opaque to KKRunner.

Declare all input, output and cache paths. Collision checks use that manifest,
including stdin, local transfer destinations and remote filesystem identity.
Sharing a calculation directory does not itself reserve the whole directory.
A missing required output is an artifact error; an absent optional output may
be skipped. Transfer errors and processing errors have distinct retry paths.

Paths written into generated SPR-KKR input are relative to the run's working
directory, even if the input file itself is elsewhere. An absolute native
``TAUDIR`` is mapped relative to the effective ``TAUPATH``. Internal open-file
paths may remain absolute to preserve ownership and cleanup guarantees.

Tracked calculations require persistent working/output paths. In particular,
``directory=False`` and temporary stderr requested by ``error_file=False`` are
rejected for tracked runs. ``error_file=None`` stores stderr beside stdout;
``error_file=True`` merges it with stdout. Temporary files for untracked runs
retain their ordinary cleanup behavior.

Do not delete files after merely losing observation. Required application
processing can resume after reconnect and needs the captured artifacts.
Processing callbacks should be idempotent or transactional because an explicit
retry after a crash can repeat side effects; exactly-once side effects are not
promised.

Remote commands for users
---------------------------

.. code-block:: console

   ase2sprkkr remote status HANDLE
   ase2sprkkr remote follow HANDLE
   ase2sprkkr remote cancel HANDLE

These commands act on the host where they run and delegate to KKRunner.
They currently accept a target-side execution handle, not a project run ID.
Internal submit/recovery/acknowledgement and MPI commands belong to
``kkrunner``. See the
`KKRunner command-line guide <https://ase2sprkkr.github.io/kkrunner/command-line.html>`_
and `recovery guide <https://ase2sprkkr.github.io/kkrunner/persistent-runs.html>`_
for the generic protocol and lifecycle contract.

Architecture for developers
-----------------------------

``sprkkr._preparation`` resolves application state and materializes files.
``PreparedCalculation`` owns prepared resources and adapts them to an execution
specification. KKRunner runs the command and persists generic execution state.
``outputs`` resolves task-owned readers and constructs ``TaskResult`` objects.
There is no scientific parser registry inside KKRunner.

KKRunner's packages are ``runs``, ``backends``, ``remote``, ``projects``, ``files``
and ``support``. Its principal classes are importable from ``kkrunner``;
implementation-level imports use their defining subpackage modules. See the
`KKRunner architecture <https://ase2sprkkr.github.io/kkrunner/architecture.html>`_
for those responsibilities and the
`API reference <https://ase2sprkkr.github.io/kkrunner/api.html>`_ for current names.
