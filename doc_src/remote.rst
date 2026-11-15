Remote running, MPI and execution profiles
==========================================

The same SPR-KKR calculations can run locally, over SSH, or through a
job scheduler. Execution profiles select the machine, executable and
execution backend without changing the scientific task.

`KKRunner documentation <https://ase2sprkkr.github.io/kkrunner/>`_
describes the execution library. Its source is on
`GitHub <https://github.com/ase2sprkkr/kkrunner>`_.

Choosing an execution profile
-----------------------------

Inspect the current settings and named profiles::

   ase2sprkkr config --show
   ase2sprkkr config --list-profiles

Without an explicitly selected named profile, ASE's ``[sprkkr]``
configuration supplies the defaults. See :doc:`installation` for
configuring local SPR-KKR executables.

For example, create a named profile for SSH execution without
installing KKRunner on the remote host (replace the paths and host)::

   ase2sprkkr config --set profiles.cluster.executor ssh
   ase2sprkkr config --set profiles.cluster.host user@cluster
   ase2sprkkr config --set profiles.cluster.backend direct
   ase2sprkkr config --set profiles.cluster.local_root /home/user/work
   ase2sprkkr config --set profiles.cluster.remote_root /scratch/user/work

Run with that profile::

   ase2sprkkr run calculation.pot --profile cluster

For an SSH job managed by a remote KKRunner supervisor, set
``profiles.cluster.backend`` to ``kkrunner`` instead. This is also
the default backend for SSH profiles that do not specify a backend.
To make the named profile the default, run
``ase2sprkkr config --set default_profile cluster``.

Select a profile from Python with
``calculator.calculate(profile="cluster")``, or construct one directly:

.. code-block:: python

   from ase2sprkkr import SPRKKR
   from ase2sprkkr.sprkkr.profile import SPRKKRProfile

   profile = SPRKKRProfile(
       executor="ssh",
       host="user@cluster",
       local_root="/path/to/local/work",
       remote_root="/path/to/remote/work",
       backend="direct",
   )
   result = SPRKKR(potential="calculation.pot").calculate(profile=profile)

Check that SSH login works, the remote directory is writable and the
remote host can run the required SPR-KKR and MPI executables. Set the
profile's ``executable_dir`` if they are not on the remote ``PATH``.

Execution backends and remote requirements
------------------------------------------

* **Direct** (``backend="direct"``): runs SPR-KKR over SSH using Python
  on the target, but **does not require KKRunner there**. It has no
  independent supervisor retaining a definitive exit status.
* **KKRunner supervisor** (``backend="kkrunner"``): runs SPR-KKR
  under a supervisor on the execution host. This requires a compatible
  KKRunner installation in the target Python environment and retains
  completion information for reconnection and recovery.
* **PSI/J** (``backend="psij"``): submits through a job scheduler,
  such as Slurm. Local scheduler access needs ``ase2sprkkr[psij]``
  on the client; scheduler access over SSH requires a compatible
  ``kkrunner[psij]`` installation on the remote host. Set
  ``psij_backend="slurm"`` for Slurm.

These are backend selections, independent of whether the calculation
is recorded in a project. See the
`KKRunner backend guide <https://ase2sprkkr.github.io/kkrunner/backends.html>`_
for further details.

Persistence, reattachment and recovery
--------------------------------------

These describe different properties of a run:

* **Persistence:** saving a run ID, its status, file paths and other
  metadata in a run store. Any backend can be used with a persistent
  project. Untracked execution is also possible with any backend.
* **Reattachment:** reconnecting to a known run after the original
  client session ends. With a direct backend, a saved process identity
  and output paths allow limited status checks and resuming observation
  of output. Once the original connection is lost, however, the final
  exit code may be unavailable, so completion cannot necessarily be
  verified as successful. The ``kkrunner`` supervisor retains
  completion information and supports fuller reattachment.
* **Submission recovery:** finding a job when the client disconnected
  before its backend ID was saved. Direct execution cannot do this.
  The supervisor or a scheduler may support it, subject to available
  receipts or scheduler records.

A lost connection does not, by itself, prove that a running process
has terminated. Equally, a persistent record alone does not guarantee
reliable reattachment or successful postprocessing. See
:doc:`execution` for enabling project tracking and loading history.

MPI and executable selection
----------------------------

MPI defaults to ``auto``. Request four ranks with:

.. code-block:: python

   from ase2sprkkr import SPRKKR

   result = SPRKKR(potential="calculation.pot", mpi=4).calculate()

A required MPI calculation must not silently fall back to serial.
If only serial executables are installed and automatic selection
cannot use them, disable MPI::

   ase2sprkkr config --set profiles.default.mpi false

For Slurm, select ``backend="psij"`` and
``psij_backend="slurm"``. With ``executor="psij"``, PSI/J
connects to a local scheduler; with ``executor="ssh"``, the
backend runs on a scheduler-accessible host via SSH.

Monitoring recorded runs
------------------------

Use ``ase2sprkkr remote status``, ``follow`` and ``cancel``
with saved run IDs. The ``--remote-id`` alternative accepts an
explicit backend handle instead. See :doc:`cli` for syntax.
