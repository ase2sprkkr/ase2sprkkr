"""
sprkkr.input_parameters.definitions - definitions of SPR-KKR tasks and their
input-parameter file schemas.

Each contained module must expose a ``task_definition`` attribute or factory
returning :class:`ase2sprkkr.input_parameters.task_definition.TaskDefinition`.
Modules are named according to the contained task (e.g. DOS, ARPES, SCF, ...).
The task definition owns execution metadata while its
``input_parameters_definition`` describes only the input-file format.
"""

import platformdirs
import os

__path__.append(os.path.join(platformdirs.user_config_dir("ase2sprkkr", "ase2sprkkr"), "input_parameters"))
