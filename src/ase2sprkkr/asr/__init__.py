""" This module just import ASR recipes.
Thus, it allows to call them in the same way as the original ASR::

  python3 -m ase2sprkkr.asr.relax --calculator "{'name':'sprkkr'}"

However, the SPRKKR calculator is available for them.
"""
try:
    import asr
    __path__ = asr.__path__
    del asr
except ImportError:
    import warnings
    warnings.warn("I am not able to import asr module, asr recipes will not be available",
                  category=RuntimeWarning)
    del warnings
