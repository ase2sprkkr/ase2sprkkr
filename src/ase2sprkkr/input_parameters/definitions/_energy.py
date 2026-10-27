""" Energy related default functions """

import functools
import numpy as np

from ..input_parameters_definitions import InputValueDefinition as V
from ...common.warnings import DataValidityError, DataValidityWarning


def validate_dos_energy(p, _values, why):
    """Validate the shared DOS/TDOBS energy-grid requirements."""
    if why == "set":
        return None
    issues = []
    for grid, ne in zip(p.ENERGY.GRID(), p.ENERGY.NE()):
        if grid in ("2", "3") and ne < 2:
            issues.append(DataValidityError("DOS GRID=2/3 requires at least two energy points."))
    if not np.isfinite(p.ENERGY.ImE()) or p.ENERGY.ImE() < 0:
        issues.append(DataValidityError("DOS ImE must be finite and nonnegative."))
    if "3" in p.ENERGY.GRID() and p.ENERGY.ImE() == 0:
        issues.append(DataValidityWarning("DOS ImE=0 puts GRID=3 on the real axis; check poles and convergence."))
    for name in ("EMIN", "EMAX", "EMINEV", "EMAXEV"):
        value = p.ENERGY[name]()
        if value is not None and not np.isfinite(value):
            issues.append(DataValidityError(f"DOS {name} must be finite."))
    return issues or None

def _energy_default(value, option):
    """Evaluate a task-supplied default, which may depend on the current mode."""
    return value(option) if callable(value) else value

def _energy_absolute_default(value, option):
    #Suppress an absolute fallback when the relative representation is active."""
    return None if option._container[option.name + 'EV']() is not None else _energy_default(value, option)


def energy_absolute_default(name, val):
    out = functools.partial(_energy_absolute_default, val)
    if callable(val):
        val = val.__doc__
    elif val is not None:
        val = str(val)
    if val is not None:
        out.__doc__ = f"{val} (if {name}_EV isn't set)"
    else:
        out.__doc__ = None

    return out


def _energy_relative_default(value, option):
    #Suppress a relative fallback when an absolute value was explicitly set.
    #Check storage with is_set(), not the absolute getter: that getter already
    #depends on this relative value, and calling it would recurse.
    return None if option._container[option.name[:-2]].is_set() else _energy_default(value, option)


def energy_relative_default(name, val):
    out = functools.partial(_energy_relative_default, val)
    if callable(val):
        val = val.__doc__
    elif val is not None:
        val = str(val)
    if val is not None:
        out.__doc__ = f"{val} (if {name} isn't set)"
    else:
        out.__doc__ = None

    return out

def energy_bound(name, values):
    """Define a paired absolute/relative bound with task-specific defaults.

    values contains (absolute Ry default, description, relative eV default).
    Either default may be a top-level callable accepting its Option. Shared
    precedence and validation are independent of the task's numerical values.
    Both absolute bounds are optional in the native input: when omitted,
    SPRKKR uses its task-specific or Fermi-energy fallback.
    """
    absolute = V(name, float, energy_absolute_default(name, values[0]) , info=values[1],
                 is_optional=True)
    relative = V(name + 'EV', float, energy_relative_default(name, values[2]),
                 info=f'{name}, given in eV with respect to the Fermi level', is_optional=True)
    return [absolute, relative]
