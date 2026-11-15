Installation and upgrading
==========================

ASE2SPRKKR provides a Python-based interface for preparing, executing and
analyzing calculations with the **separately installed** SPR-KKR software.
For the first calculation, continue with :doc:`quickstart` after installing
the interface and configuring the executables.

Requirements
------------

* **Python 3.9 or newer** for ASE2SPRKKR 4.x. Earlier releases have
  different Python requirements; version 3.6 was the last supporting Python 3.8.
* Native SPR-KKR executables, libraries and any required MPI runtime for
  **running** calculations. The interface can be used to inspect task
  definitions and parse existing results without these executables.

Get the native binaries from the
`SPR-KKR download page <https://sprkkr.org/download/>`_ and follow the
installation instructions distributed with the program. The
`SPR-KKR manual <https://sprkkr.org/manual/>`_ covers the native
build requirements, methodology and numerical settings. ASE2SPRKKR does
not distribute or compile SPR-KKR.

Check whether native executables are on ``PATH`` (for example,
``command -v kkrscf`` on Unix-like systems). The command
``ase2sprkkr info --executables`` lists known SPR-KKR programs,
not necessarily the binaries installed on your system.

Install or upgrade ASE2SPRKKR
-----------------------------

Install the latest stable release into the Python environment where
you intend to run ASE2SPRKKR::

   python -m pip install --upgrade ase2sprkkr
   ase2sprkkr --version

Pip does not select beta releases by default. To install the newest
available pre-release (for example in the 4.0 beta series), use::

   python -m pip install --upgrade --pre ase2sprkkr

The version 4 features described in this documentation may not be
available in earlier stable releases. Check the installed version with
``ase2sprkkr --version`` and consult :doc:`changelog`.

The Python distribution automatically installs its required packages,
including KKRunner on the local machine. To reproduce a specific
version, pin it with
``python -m pip install 'ase2sprkkr==VERSION'``.

Use ``python -m pip`` with the intended Python interpreter. Do not use
``sudo pip`` on distribution-managed Python installations. If your
system restricts installation into the main environment, use the
:ref:`installation-venv` below.

Configure SPR-KKR executables
-----------------------------

If SPR-KKR executables are already on ``PATH``, no additional setting is
needed. Otherwise, specify the directory containing the executables::

   ase2sprkkr config --set profiles.default.executable_dir /opt/sprkkr/bin

Replace the example path with the real directory. The command edits the
active execution profile; by default this is ASE's ``[sprkkr]``
configuration. To edit ASE's profile explicitly, use ``profiles.ase``.

If your executables use a suffix (for example ``kkrscf_v8``), configure
it as well::

   ase2sprkkr config --set profiles.default.executable_suffix _v8

Inspect the current configuration and existing profiles with::

   ase2sprkkr config --show
   ase2sprkkr config --list-profiles

The default MPI policy is ``auto``. When only serial executables are
installed and automatic selection cannot use them, explicitly disable
MPI with ``ase2sprkkr config --set profiles.default.mpi false``.
Execution profiles, MPI and clusters are covered in :doc:`remote`.

Named profiles use ASE2SPRKKR's YAML configuration; the ASE profile
uses ASE's own ``config.ini``. CLI ``config --set`` changes are
saved immediately. Changes made through Python
``config.profiles.ase.OPTION.set(value)`` stay in memory until saved
with ``config.save_user_preferences()``.

.. _installation-venv:

Virtual environments
--------------------

A virtual environment keeps an ASE2SPRKKR installation and its Python
dependencies separate. This is useful when testing a different version,
avoiding package conflicts, or working on a system-managed Python
installation. It is optional if your current environment is suitable.

For example, on Linux/macOS::

   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade ase2sprkkr
   ase2sprkkr --version

On Windows, activate with ``.venv\Scripts\activate``.
Run ``deactivate`` to leave the environment. Installing the Python
interface in a venv does not install the native SPR-KKR programs;
configure or find the same native executables as described above.

.. _installation-upgrade:

Upgrading from 3.x or older versions
------------------------------------

Upgrade the package in your regular Python environment with::

   python -m pip install --upgrade ase2sprkkr

When upgrading an installation with legacy Python configuration
(``~/.config/ase2sprkkr/__init__.py`` on typical Linux systems), run::

   ase2sprkkr migrate

This command converts supported declarative configuration settings
to ``config.yaml``, including legacy execution settings migrated into
a named execution profile. It leaves unrelated Python statements alone.
By default it removes the migrated assignments from the old hook;
use ``ase2sprkkr migrate --keep-init`` to retain that hook unchanged,
or ``--delete-init`` to remove it. The existing
``ase2sprkkr config --migrate`` spelling is also supported.
If the destination YAML file already exists, the migration does not
overwrite it. Review the settings afterward with
``ase2sprkkr config --show``.

For an **isolated trial upgrade**, follow :ref:`installation-venv`
instead of changing the active environment.

Version 4 changes task definitions, execution and result handling.
When upgrading existing scripts, check the following:

* Task names and options: ``ase2sprkkr info --task`` and :doc:`usage`.
* Execution profiles and SSH/schedulers: :doc:`remote`.
* Persistent projects and run history: :doc:`execution`.
* Result parsing, plotting and NeXus export: :doc:`result-files`.

See :doc:`changelog` for the complete release-specific changes.
Test important calculations before migrating production workflows.

.. _installation-optional:

Optional dependencies
---------------------

The core Python installation does not need all optional plotting, data-export
or integration dependencies. Add only the features you use; extras
can be combined, for example::

   python -m pip install --upgrade 'ase2sprkkr[plot,nexus]'

The extras are:

``plot``
   Matplotlib, SciPy and scikit-image for supported scientific plots
   and Grace ``.agr`` files. The declared extra is currently called
   ``plot``; use the name listed by your installed version.

``nexus``
   ``h5py`` and ``nexusformat`` for ARPES SPC to NeXus/HDF5
   export. See :doc:`result-files`.

``mendeleev``
   Extended element metadata, accessed lazily through
   ``AtomicType.mendeleev``. Basic atomic information uses ASE or
   built-in tables without requiring this extra.

``psij``
   ``kkrunner[psij]`` for integration with PSI/J job schedulers.
   Remote scheduler requirements are explained in :doc:`remote`.

``jupyter``
   JupyterLab, nbformat and nbconvert for notebook workflows.

``es_finder``
   Integration with the **proprietary, non-public** ``es_finder``
   empty-sphere finder. This package is not distributed through PyPI:
   ``ase2sprkkr[es_finder]`` installs only the public
   ``pymatgen`` prerequisite, **not the proprietary program**.
   The integration requires separately obtaining and installing
   the proprietary software.

The `ASE2SPRKKR GUI <https://github.com/ase2sprkkr/ase2sprkkr-gui>`_
is a separately installed application, not an ``ase2sprkkr[gui]``
extra. Install optional packages in the Python environment used by
the corresponding application.

Developer installation
----------------------

The repository uses separate branches for ongoing work and releases:

* ``develop`` — active development and the normal base for contributions.
  It may contain code not yet included in a published version.
* ``master`` — the beta-release line; it may lag behind ``develop``.
* ``release`` — the stable-release line, when present; it tracks
  published stable versions rather than ongoing development.

Use a version **tag** for a reproducible release checkout instead of
depending on the moving branch heads. For development, clone
``develop`` and install in editable mode::

   git clone https://github.com/ase2sprkkr/ase2sprkkr.git
   cd ase2sprkkr
   git switch develop
   python -m pip install -e .

Source builds require a native compiler and CMake. Pip's isolated
build installs the declared Python build dependencies (including
scikit-build-core and Cython); using ``--no-build-isolation`` means
installing them yourself. For a smoke test without native executables::

   ase2sprkkr test --no-kkr --fast

See :doc:`philosophy` for developer architecture and :doc:`remote`
for SSH and scheduler requirements.

Building this documentation
---------------------------

The Sphinx sources live in ``doc_src/``; generated HTML lives in
``docs/``. The authoritative GitHub README is
``doc_src/README.rst``, synchronized to ``README.md`` by
``make docs`` or ``make doc-readme``. Verify the synchronization
using ``make doc-readme-check``.

From an editable source checkout, install the documentation tools::

   python -m pip install sphinx sphinx-toolbox sphinx-rtd-theme matplotlib

``make doc-build`` generates the HTML and source-driven example gallery.
``make docs`` additionally regenerates ``README.md`` from
``doc_src/README.rst``; this requires the external ``pandoc`` program.
To verify that the checked-in Markdown matches its source, run
``make doc-readme-check`` (also requiring Pandoc).

The full HTML build imports ASE2SPRKKR and its dependencies; install
the package from this checkout before invoking Sphinx.
