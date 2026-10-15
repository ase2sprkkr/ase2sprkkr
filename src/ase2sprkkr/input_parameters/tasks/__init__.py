"""Explicit task definitions for tasks sharing input schemas.

Only non-1:1 task/schema relationships belong here. Standard tasks are exposed
through their InputParametersDefinition._task_definition.
"""

import os
import platformdirs

__path__.append(
    os.path.join(
        platformdirs.user_config_dir("ase2sprkkr", "ase2sprkkr"),
        "tasks",
    )
)
