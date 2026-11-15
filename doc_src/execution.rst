Persistent project calculations
===============================

Calculations are untracked by default. A persistent project stores run
records for tracking progress and loading results later. This is
independent of the execution backend. Reconnecting to a known running
job is a separate capability: direct execution supports limited output
reattachment without guaranteeing the lost exit code, whereas the
``kkrunner`` supervisor retains completion information. See
:doc:`quickstart` for ordinary calculations and :doc:`remote`
for the differences between persistence, reattachment and recovery.

Create and track a run
----------------------

Create a project once and pass its ID to the calculation:

.. code-block:: python

   from ase.build import bulk
   from ase2sprkkr import SPRKKR
   from ase2sprkkr.sprkkr.projects import default_project_repository

   project = default_project_repository().create("Aluminium")
   calculator = SPRKKR(atoms=bulk("Al", "fcc", a=4.05), directory="al-scf")
   run = calculator.create_run(project=project.id)
   print("Project ID:", project.id)
   print("Run ID:", run.id)
   run.submit()
   result = run.wait()

``create_run()`` prepares the calculation, ``submit()`` starts execution,
and ``wait()`` returns the SPR-KKR result. The project ID identifies
the project; the run ID identifies the calculation. The output files
must remain available for loading and postprocessing.

For a synchronous calculation with project tracking, use
``calculator.calculate(project=project.id)`` instead of
``create_run()``, ``submit()`` and ``wait()``.

On the command line, select an existing project by UUID and optionally
an execution profile::

   ase2sprkkr run example.pot --task SCF --project-id PROJECT_UUID --profile cluster

The command prints the run ID after submission and normally waits for the
result. To submit a project calculation without waiting, use::

   ase2sprkkr run example.pot --task SCF --project-id PROJECT_UUID --profile cluster --submit-only

``--submit-only`` requires ``--project-id`` and a managed ``kkrunner`` or
``psij`` backend. A tracked direct run must remain under CLI observation
until completion. The working directory defaults to the project's directory
when configured, or to the current directory. Use ``--directory PATH`` to
override it. See :doc:`cli` for all options.

Reconnect or load history
-------------------------

After restarting the client, reopen the project manager to resume
observing a known nonterminal run, where backend capabilities permit:

.. code-block:: python

   from ase2sprkkr.sprkkr.projects import project_execution_manager
   from ase2sprkkr.outputs.project_results import attach_project_run

   manager = project_execution_manager(project="PROJECT_UUID")
   run = attach_project_run(manager, "RUN_ID")
   result = run.wait()

For a completed historical run, load the parsed result rather than
reattaching the process:

.. code-block:: python

   from ase2sprkkr.outputs.task_result import TaskResult

   result = TaskResult.from_run("RUN_ID", manager=manager)

Replace the placeholder IDs with those printed when the project and run
were created. History requires access to the same project run store and
retained output files. Failed or interrupted postprocessing may be
retried without resubmitting SPR-KKR using
``retry_project_postprocessing(manager, run_id)``.

Inspect and follow a run
------------------------

.. code-block:: console

   ase2sprkkr remote status RUN_ID
   ase2sprkkr remote follow RUN_ID --info --from-start
   ase2sprkkr remote cancel RUN_ID

These commands use the saved project run ID. See :doc:`cli`
for other options, or :doc:`remote` for backend requirements.
