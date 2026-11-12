CLI command plugins
===================

An installed Python distribution can extend the ``ase2sprkkr`` command
without modifying ASE2SPRKKR. A plugin exposes a **Python package** through
the entry point named ``cli`` in the group ``ase2sprkkr``:

.. code-block:: toml

   [project.entry-points."ase2sprkkr"]
   cli = "my_plugin.commands"

The target must be an importable package (with an ``__init__.py``), not a
single command module. Its immediate Python modules become CLI commands:

.. code-block:: text

   my_plugin/
       __init__.py
       commands/
           __init__.py
           analyse.py
           show_results.py

For example, ``show_results.py`` makes ``ase2sprkkr show-results``
available. Module filenames use underscores; CLI names use hyphens. Modules
without a ``parser`` function are ignored.

Each command has the **same interface as native commands**:

.. code-block:: python

   help = "Analyse a calculation."
   description = "Detailed description of the plugin command."
   aliases = ["inspect"]  # optional

   def parser(parser):
       parser.add_argument("filename")

   def run(args, global_args):
       from my_plugin.analysis import analyse  # expensive imports belong here

       analyse(args.filename)

The required module members are ``help``, ``description``,
``parser(parser)`` and ``run(args, global_args)``. Optional members
``aliases``, ``extra_help`` and ``unknowns`` behave as they do for
native commands. ``global_args`` contains the global ``debug`` and
``profile`` flags. The plugin's argument parser is initialized when its
command (or alias) is selected.

Install the plugin distribution in the same Python environment as
ASE2SPRKKR, for example with ``pip install my-plugin``. Its commands
then appear in ``ase2sprkkr --help``.

Native command names and aliases take precedence over plugin commands.
If a plugin claims an already registered name or alias, that **command**
is skipped with a diagnostic on stderr; other plugin commands remain
available. Plugins are processed in deterministic entry-point target
order. A plugin whose entry point or command module fails to import is
also reported and skipped. Command modules should avoid expensive
top-level imports so that normal CLI startup remains lightweight.
