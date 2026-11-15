ASE2SPRKKR architecture for developers
======================================

ASE2SPRKKR separates scientific definitions and result parsing from
the execution of native SPR-KKR programs. The ``SPRKKR`` calculator
is the main user-facing interface; scientific file formats and
execution profiles are also available independently.

Configuration files and definitions
-----------------------------------

SPR-KKR input and potential files use configuration objects.
Instances hold sections and option values; definitions specify
grammar, types, defaults, validation and parsing/writing rules.
``ase2sprkkr.common`` provides ``ConfigurationContainer``,
``Option``, ``ContainerDefinition`` and ``ValueDefinition``,
together with grammar types and
`pyparsing <https://pyparsing-docs.readthedocs.io/en/latest/>`_.

TaskDefinition: task identity versus input schema
-------------------------------------------------

``ase2sprkkr.input_parameters.task_definition.TaskDefinition``
is an **immutable description of one SPR-KKR task**. It declares
the canonical ``name``, associated ``input_parameters_definition_name``,
native ``executable``, MPI support, output reader, experimental status
and ``TaskFile`` input/output/cache manifests. It is not a second
parser or a holder of the parameter values.

The distinction matters because multiple tasks can share one input
parameter schema. For example, the DMI task uses JXC input parameters:

.. code-block:: python

   from ase2sprkkr.input_parameters.task_definition import TaskDefinition

   def task_definition():
       return TaskDefinition(
           "JXC", name="DMI", executable="kkrgen",
           mpi=True, experimental=True, output_reader="JXC",
       )

A task with its own input schema can declare
``schema._task_definition`` alongside the corresponding
``input_parameters_definition()`` function, as DOS and SCF do.
Tasks sharing an existing schema are instead defined separately
in ``ase2sprkkr.input_parameters.tasks``.

**Resolve task identities through the TaskRegistry**, not by assuming
that an input schema name identifies exactly one task.
``InputParametersRegistry`` owns input schema resolution.
The public API ``InputParameters.for_task(task_name)`` creates
task-bound input parameters, and ``ase2sprkkr info --task``
shows the registered tasks.

When adding a task, define its native executable and supported
``TaskFile`` manifests, select its parameter schema and output
reader, and add focused tests. Avoid maintaining parallel task
identities or hard-coding known files in the calculator.

Scientific layers
-----------------

* ``ase2sprkkr.input_parameters`` — task and schema registries,
  task definitions and native input parsing.
* ``ase2sprkkr.potentials`` — potential-file representations.
* ``ase2sprkkr.sprkkr`` — ASE calculator, atom/site properties,
  preparation, profiles and symmetry.
* ``ase2sprkkr.outputs`` — task output readers, ``TaskResult``.
* ``ase2sprkkr.result_files`` — structured result readers, Grace
  ``.agr`` and NeXus export.
* ``ase2sprkkr.bindings`` — optional external integrations.
* ``ase2sprkkr_cli`` — lightweight command dispatcher and plugins.

Execution boundary
------------------

Preparation via ``sprkkr._preparation`` materializes native
inputs and file manifests. ``PreparedCalculation`` adapts those
to KKRunner specifications.

`KKRunner <https://ase2sprkkr.github.io/kkrunner/>`_ handles
local and MPI execution, SSH, schedulers, transfers, attachment
and persistent runs. It does not depend on ASE2SPRKKR's physics
parsers. ASE2SPRKKR interprets completed calculations through
its output reader layer. See :doc:`remote` for execution profiles and :doc:`execution`
for persistent calculations.

Development and documentation
-----------------------------

Source builds and optional tools are described in :doc:`installation`.
The authoritative README is ``doc_src/README.rst`` and
``make docs`` synchronizes the root ``README.md``. API documentation
and the examples gallery are generated from Python sources; avoid
editing generated pages. See :doc:`cli-plugins` for CLI extension.

Files and processing (maintainers)
----------------------------------

Task definitions provide KKRunner with input, output and cache file
declarations. Scientific format descriptions and result readers belong to
ASE2SPRKKR, not to the execution library. The files needed for
reconnection must remain available until postprocessing finishes.
A required output must be present; optional outputs can be absent.
Application postprocessing may be retried after interruption, so its
callbacks should be safe to repeat. See the
`KKRunner architecture <https://ase2sprkkr.github.io/kkrunner/architecture.html>`_
before changing this execution boundary.
