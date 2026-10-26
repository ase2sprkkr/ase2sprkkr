""" Energy related default functions """

import functools

from ..input_parameters_definitions import InputValueDefinition as V

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
