""" This module just import ASR recipes.
Thus, it allows to call them in the same way as the original ASR::

  python3 -m ase2sprkkr.asr.relax --calculator "{'name':'sprkkr'}"

However, the SPRKKR calculator is available for them.
"""

import asr
__path__ = asr.__path__
del asr
