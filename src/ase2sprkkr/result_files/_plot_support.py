"""Lazy delegates keeping optional plotting libraries out of parsers."""
from functools import partial, update_wrapper


class Multiplot:
    """Construct the real plotting context only when explicitly requested."""

    def __new__(cls, *args, **kwargs):
        from ase2sprkkr.gui.plot import Multiplot as RealMultiplot
        return RealMultiplot(*args, **kwargs)


def change_default_kwargs(func, **kwargs):
    """Bind plotting defaults without importing Matplotlib."""
    wrapper = partial(func, **kwargs)
    update_wrapper(wrapper, func)
    return wrapper


def single_plot(*args, **kwargs):
    from ase2sprkkr.gui.plot import single_plot as real
    return real(*args, **kwargs)


def colormesh(*args, **kwargs):
    from ase2sprkkr.gui.plot import colormesh as real
    return real(*args, **kwargs)


def lineplot(*args, **kwargs):
    from ase2sprkkr.gui.plot import lineplot as real
    return real(*args, **kwargs)


def set_up_common_plot(*args, **kwargs):
    from ase2sprkkr.gui.plot import set_up_common_plot as real
    return real(*args, **kwargs)
