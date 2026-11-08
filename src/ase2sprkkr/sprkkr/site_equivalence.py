"""Define the exact equivalence classes of sites passed to spglib."""

import numpy as np

from ..common.unique_values import UniqueValuesMapping
from .occupations import Occupation


def occupation_key(occupation):
    """Exact chemical composition, including vacancies and independent of order."""
    return tuple(sorted(
        (symbol, float(amount)) for symbol, amount in occupation.to_dict(True).items()
    ))


def site_equivalence_tags(
    atoms, *, equivalence="site-types", atomic_kinds=None,
    assigned_sites=None, region_memberships=None,
):
    """Classify sites once, before the symmetry search.

    Uninitialized sites use ASE occupations. Initialized sites in strict mode
    retain their SiteType identities and independent per-site magnetisations.
    Explicit atomic_kinds override the species classification.
    assigned_sites preserves SiteType sharing established by earlier region
    calculations. region_memberships distinguishes atoms in different regions
    even if they otherwise have identical SiteTypes.
    """
    if equivalence not in ("site-types", "occupation"):
        raise ValueError("equivalence must be 'site-types' or 'occupation'")

    arrays = atoms.arrays
    sites = arrays["sprkkr_sites_data"] if "sprkkr_sites_data" in arrays else None
    kinds = arrays["spacegroup_kinds"] if "spacegroup_kinds" in arrays else None
    initial_magmoms = (
        arrays["initial_magmoms"] if "initial_magmoms" in arrays else None
    )
    if atomic_kinds is not None and len(atomic_kinds) != len(atoms):
        raise ValueError("atomic_kinds must contain one value per site")
    if assigned_sites is not None and len(assigned_sites) != len(atoms):
        raise ValueError("assigned_sites must contain one value per site")
    if region_memberships is not None and len(region_memberships) != len(atoms):
        raise ValueError("region_memberships must contain one value per site")

    partial = sites is not None and any(not s for s in sites)
    keys = []
    for i in range(len(atoms)):
        site = sites[i] if sites is not None and sites[i] else None
        if atomic_kinds is not None:
            # Explicit atomic kinds completely replace implicit site classes.
            kind = atomic_kinds[i]
            key = tuple(kind.tolist()) if isinstance(kind, np.ndarray) else kind
        elif site is not None and equivalence == "site-types" and not partial:
            key = ("site-type", site.site_type, tuple(site.magnetisation))
        else:
            occ = site.occupation if site is not None else Occupation.for_atom_from_ase_atoms(atoms, i)
            key = ("occupation", occupation_key(occ))

        if atomic_kinds is None and equivalence == "site-types":
            if kinds is not None:
                key = (key, int(kinds[i]))
            if initial_magmoms is not None:
                moment = np.asarray(initial_magmoms[i]).reshape(-1)
                key = (key, tuple(float(component) for component in moment))
        if assigned_sites is not None:
            # Group by the shared SiteType, not by each distinct Site object.
            site = assigned_sites[i]
            key = (key, site.site_type if site else None)
        if region_memberships is not None:
            key = (key, tuple(region_memberships[i]))
        keys.append(key)

    return UniqueValuesMapping.from_values(keys).mapping.astype(np.intc, copy=False)


def refine_partial_site_tags(tags, dataset, sites):
    """Refine chemical symmetry using known SiteTypes of partially set sites.

    Uninitialized sites inherit the known identity of their preliminary
    symmetry orbit. Distinct known SiteTypes or site magnetisations split
    an orbit; no Site objects are created or changed.
    """
    known = {}
    for i, site in enumerate(sites):
        if site:
            known.setdefault(
                int(dataset.equivalent_atoms[i]),
                (site.site_type, tuple(site.magnetisation)),
            )

    keys = [
        (
            int(tag),
            (site.site_type, tuple(site.magnetisation))
            if site else known.get(int(dataset.equivalent_atoms[i])),
        )
        for i, (tag, site) in enumerate(zip(tags, sites))
    ]
    refined = UniqueValuesMapping.from_values(keys).mapping
    if UniqueValuesMapping.are_equivalent(tags, refined):
        return tags
    return refined
