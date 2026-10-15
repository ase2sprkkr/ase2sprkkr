"""Input-parameter schema definitions.

Each public module exposes ``input_parameters_definition``. Standard schemas
attach their 1:1 task as ``_task_definition``; schemas shared by several tasks
leave it unset and resolve task names from concrete parameter values.
"""

import os
import platformdirs

__path__.append(
    os.path.join(
        platformdirs.user_config_dir("ase2sprkkr", "ase2sprkkr"),
        "input_parameters",
    )
)
