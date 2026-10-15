from ase2sprkkr.sprkkr.sites import Site, SiteType
from ase2sprkkr.sprkkr.sprkkr_atoms import SPRKKRAtoms


def test_site_type_break_symmetry_detaches_all_shared_sites():
    atoms = SPRKKRAtoms("Fe")
    shared = SiteType(atoms, {"Fe": 1.0})
    sites = [Site(shared) for _ in range(3)]

    assert shared.sites == set(sites)
    assert all(site.has_symmetry for site in sites)

    shared.break_symmetry()

    site_types = [site.site_type for site in sites]
    assert len({id(site_type) for site_type in site_types}) == len(sites)
    assert all(not site.has_symmetry for site in sites)
