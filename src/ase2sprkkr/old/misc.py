# coding: utf-8
# vim: set et sw=4 ts=4 sts=4 nu ai cc=80 fdm=indent mouse=a:

"""
Module misc
===========

This module is a place to put everything that could not fit
anywhere else...

"""

import logging
import numpy as np

from ase.spacegroup import get_spacegroup


LOGGER = logging.getLogger('sprkkr')


def set_log_level(level):
    lvl = getattr(logging, level.upper())
    LOGGER.setLevel(lvl)


def set_occupancy(atoms, site, symbol, concentration):
    sites = atoms.get_scaled_positions()
    n = len(sites)
    if site >= n:
        raise ValueError(f"The site number should be <= {n:d}")
    # get the previously defined occupancy if it exists or create it
    occupancy = atoms.info.get('occupancy', {})
    site_occ = occupancy.get(site, {})
    site_occ.update({symbol: concentration})
    occupancy.update({site: site_occ})
    atoms.info.update({'occupancy': occupancy})


def get_occupancy(atoms):
    occupancy = atoms.info.get('occupancy', {})
    sites = atoms.get_scaled_positions()
    for site in range(len(sites)):
        site_occ = occupancy.get(site, {})
        if not site_occ:
            # find coords of site
            site_xyz = sites[site]
            # get the atom at such coordinates
            i = np.where((atoms.get_scaled_positions() == sites[site]).all(axis=1))[0][0]
            atom = atoms[i]
            site_occ = {atom.symbol: 1.}
        occupancy.update({site: site_occ})
    return occupancy


class AttrDict(dict):
    """
    A dict with the attribute access to its items.
    """

    def __getattr__(self, name):
        try:
            return self[name]

        except KeyError:
            raise AttributeError(name)

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __str__(self):
        if self.keys():
            return '\n'.join(
                [self.__class__.__name__ + ':'] +
                self.format_items()
            )

        else:
            return self.__class__.__name__

    def __repr__(self):
        return self.__class__.__name__

    def __dir__(self):
        return list(self.keys())

    def copy(self):
        return type(self)(self)

    def format_items(self):
        num = max(map(len, list(self.keys()))) + 1
        return [key.rjust(num) + ': ' + repr(val)
                for key, val in sorted(self.items())]
