"""
SPR-KKR task definitions.

Each task module exports a ``task_definition()`` factory. The factory returns a
:class:`~ase2sprkkr.input_parameters.task_definition.TaskDefinition`, which
owns both task execution semantics and the corresponding input-file schema.
Modules are discovered lazily with :mod:`pkgutil`; user task modules can be
placed in the configured input-parameters directory as before.
"""

import platformdirs
import os

__path__.append(
    os.path.join(
        platformdirs.user_config_dir("ase2sprkkr", "ase2sprkkr"),
        "input_parameters",
    )
)
