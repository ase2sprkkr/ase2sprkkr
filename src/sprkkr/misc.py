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


logging.basicConfig(level=logging.DEBUG)
LOGGER = logging.getLogger('sprkkr')


def set_log_level(level):
    lvl = getattr(logging, level.upper())
    LOGGER.setLevel(lvl)


def set_occupancy(atoms, site, symbol, concentration):
    # get the spacegroup
    sg = get_spacegroup(atoms)
    # get the number of unique sites
    sites = sg.unique_sites(atoms.positions)
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
    sg = get_spacegroup(atoms)
    sites = sg.unique_sites(atoms.get_scaled_positions())
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



