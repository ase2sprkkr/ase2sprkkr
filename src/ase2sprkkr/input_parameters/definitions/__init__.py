"""SPR-KKR task definitions and their input-parameter schemas.

Bundled task modules are named after the task they define. User task modules are
also discovered from the ``input_parameters`` directory under the ase2sprkkr
user configuration directory appended to ``__path__`` below.

Each public user task module must expose a ``task_definition`` attribute or a
zero-argument factory returning
:class:`ase2sprkkr.input_parameters.task_definition.TaskDefinition`. Starting
with version 3.6 this replaces the legacy ``input_parameters`` discovery
export. Helper modules should use private names beginning with ``_``.
"""

import os

import platformdirs


__path__.append(
    os.path.join(
        platformdirs.user_config_dir("ase2sprkkr", "ase2sprkkr"),
        "input_parameters",
    )
)
