"""Temporary internal bridge for persistence code during lifecycle migration.

The public enum contains the new lifecycle only.  The legacy persistence module
still spells three historical names internally; expose them as ordinary class
attributes (not Enum members) until that module is mechanically cleaned up.
"""

from .base import FailurePhase, RunState

RunState.FINALIZING = RunState.POSTPROCESSING
RunState.FINALIZATION_FAILED = RunState.FAILED
FailurePhase.FINALIZATION = FailurePhase.POSTPROCESSING
