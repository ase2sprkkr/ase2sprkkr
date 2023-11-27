"""
sprkkr.input_parameters.definitions - definitions of InputParameterss: definition of configurtion options and their
default values for SPR-KKR input_parameters.

Each contained module should have input_parameters attribute, that holds the definition of a task. The modules
should be named according to the name of the contained task (e.g. DOS, ARPES, SCF, ...)
"""
import platformdirs
import os

__path__.append(
      os.path.join(platformdirs.user_config_dir('ase2sprkkr', 'ase2sprkkr'), 'input_parameters')
)
