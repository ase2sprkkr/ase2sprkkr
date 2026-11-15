Command-line reference
======================

The ``ase2sprkkr`` command is installed by the Python distribution.
Use the installed version's help for the authoritative options::

   ase2sprkkr --help
   ase2sprkkr --version
   ase2sprkkr <command> --help

Some features require separately installed SPR-KKR binaries or
additional Python packages. See :doc:`installation`.

Discover tasks, formats and executables
---------------------------------------

.. code-block:: console

   ase2sprkkr info --task
   ase2sprkkr info --task=SCF
   ase2sprkkr info --task=SCF.ENERGY
   ase2sprkkr info --output-file
   ase2sprkkr info --executables

``info --task`` lists the available tasks, marking experimental ones.
Use ``--verbose`` for additional parameter descriptions.

Configuration and examples
--------------------------

.. code-block:: console

   ase2sprkkr config --show
   ase2sprkkr config --path
   ase2sprkkr config --list-profiles
   ase2sprkkr config --info
   ase2sprkkr config --set profiles.default.executable_dir /opt/sprkkr/bin
   ase2sprkkr config --set profiles.default.executable_suffix _v8
   ase2sprkkr migrate
   ase2sprkkr example
   ase2sprkkr example --regex ARPES
   ase2sprkkr example 8 --copy ./my-example

Use ``profiles.default`` to configure the active profile. By default,
this is ASE's ``[sprkkr]`` configuration; ``profiles.ase`` edits it
explicitly. Named ASE2SPRKKR profiles use YAML configuration.
``ase2sprkkr migrate`` converts supported older settings.
For further profile options see :doc:`remote`; for upgrading and
configuration storage see :doc:`installation`.
``config --show-default`` displays the available YAML settings.

Run and inspect files
---------------------

.. code-block:: console

   ase2sprkkr run example.pot --task SCF
   ase2sprkkr run example.pot --task SCF --profile cluster --project-id PROJECT_UUID
   ase2sprkkr run example.pot --project-id PROJECT_UUID --profile cluster --submit-only
   ase2sprkkr run example.pot --project-id PROJECT_UUID --directory ./work
   ase2sprkkr show-structure example.pot --no-view -o example.cif
   ase2sprkkr symmetry input.cif
   ase2sprkkr symmetry input.cif --primitive -o primitive.cif
   ase2sprkkr k-path example.pot

``run`` accepts a potential file or an input file and optional
``NAME=VALUE`` overrides. ``--profile`` selects a named execution
profile; ``--project-id`` creates a persistent calculation in an existing
project and prints its run ID. Without ``--project-id`` the calculation
is untracked. The command normally waits for completion. For a managed
``kkrunner`` or ``psij`` backend, ``--submit-only`` with ``--project-id``
submits the run and exits without waiting. ``--directory`` overrides the
project's working directory, if any. See :doc:`execution` for project
tracking and :doc:`remote` for backends. The ``run`` command's
``--profile NAME`` is distinct from the top-level ``--cprofile``
option for profiling the Python command. ``show-structure`` (alias ``view``) can display or
export the structure, and ``--close-after`` can close its viewer
automatically. The interactive ``k-path`` tool requires
a graphical environment. Symmetry analysis and conversion are described
in :doc:`symmetry`.

Plotting, Grace and NeXus
-------------------------

.. code-block:: console

   ase2sprkkr plot calculation.dos
   ase2sprkkr plot calculation.dos --show-values
   ase2sprkkr plot calculation.dos -o dos.png
   ase2sprkkr plot plot.agr --show-values
   ase2sprkkr plot plot.agr -o grace.png
   ase2sprkkr nexus arpes.spc -o arpes.nxs
   ase2sprkkr nexus arpes.spc -o arpes.nxs --input-file arpes.inp

Install ``ase2sprkkr[plot]`` for plotting or
``ase2sprkkr[nexus]`` for NeXus export; see
:ref:`installation-optional`. The ``nexus`` command takes a positional ``.spc`` filename;
``-o/--output`` selects the destination and ``-i/--input-file``
supplies optional provenance.
See :doc:`result-files` for supported data and limitations.

Remote status, output and cancellation
--------------------------------------

.. code-block:: console

   ase2sprkkr remote status RUN_ID
   ase2sprkkr remote status RUN_ID --status-only
   ase2sprkkr remote follow RUN_ID --info --from-start
   ase2sprkkr remote follow RUN_ID --regex 'ERROR'
   ase2sprkkr remote cancel RUN_ID

These default to a **persistent project run ID** in the current user's
project run store. To address an explicit target-side backend handle,
use ``--remote-id`` **instead** of the positional run ID::

   ase2sprkkr remote status --remote-id HANDLE
   ase2sprkkr remote follow --remote-id HANDLE --from-start

``remote follow --raw`` emits binary framed output; do not combine
it with ``--info`` or ``--regex``. ``--stdout-path`` and
``--stderr-path`` apply only with ``--remote-id``.
See :doc:`execution` for persistence and :doc:`remote` for profiles.

Other tools
-----------

``ase2sprkkr uppasd`` can write inputs for UppASD from JXC/DMI
results; ``ase2sprkkr nomad`` provides NOMAD integration;
``ase2sprkkr shell`` opens an interactive environment. For local
sanity checks without installed SPR-KKR executables::

   ase2sprkkr test --no-kkr --fast

Third-party distributions can add commands through documented entry
points; see :doc:`cli-plugins`.
