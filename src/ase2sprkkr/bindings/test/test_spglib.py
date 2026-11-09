"""Single spglib entry point and original-basis geometry idealization."""

from types import SimpleNamespace

import numpy as np
import pytest
from ase.build import bulk

from ..spglib import (
    TransformedCell, idealize_cell, primitive_site_groups,
    symmetry_dataset, transform_cell,
)
from ...sprkkr.site_equivalence import site_equivalence_tags
from ...sprkkr.sprkkr_atoms import SPRKKRAtoms
from ...sprkkr.atoms_region import AtomsRegion


def test_dataset_and_transformation_use_tag_classes():
    atoms = bulk("NaCl", "rocksalt", a=5.64)
    tags = site_equivalence_tags(atoms)
    dataset = symmetry_dataset(atoms.cell, atoms.get_scaled_positions(), tags)
    assert dataset.number == 225
    transformed = transform_cell(
        atoms.cell, atoms.get_scaled_positions(), tags, to_primitive=True,
    )
    assert len(transformed.tags) == 2
    assert set(transformed.tags) == set(tags)


def test_site_tags_respect_initialized_site_types():
    atoms = SPRKKRAtoms(bulk("Fe", "bcc", a=2.87, cubic=True), symmetry=False)
    _ = atoms.sites
    strict = site_equivalence_tags(atoms)
    relaxed = site_equivalence_tags(atoms, equivalence="occupation")
    assert strict[0] != strict[1]
    assert relaxed[0] == relaxed[1]


def test_region_wrapped_arrays_and_partial_sites():
    atoms = SPRKKRAtoms(bulk("Fe", "bcc", a=2.87, cubic=True), symmetry=False)
    _ = atoms.sites
    region = AtomsRegion("bulk", slice(None), atoms.cell.array, [True, True, True], atoms=atoms)
    assert len(site_equivalence_tags(region)) == len(atoms)
    sites = atoms.sites.copy()
    sites[1] = 0
    atoms.set_array("sprkkr_sites_data", sites)
    assert len(site_equivalence_tags(atoms)) == len(atoms)
    assert symmetry_dataset(
        atoms.cell, atoms.get_scaled_positions(), site_equivalence_tags(atoms)
    )


def test_idealize_tetragonal_metric_preserves_distinct_c():
    rotations = [
        np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float),
        np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]], dtype=float),
        np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]], dtype=float),
        np.eye(3),
    ]
    dataset = SimpleNamespace(
        rotations=rotations, translations=np.zeros((4, 3))
    )
    lattice = np.diag([4.0, 4.002, 6.0])
    result, changed = idealize_cell(
        lattice, [[0, 0, 0]], [1], dataset, symprec=0.01
    )
    lengths = np.linalg.norm(result.cell, axis=1)
    assert changed
    assert lengths[0] == pytest.approx(lengths[1], abs=1e-10)
    assert lengths[2] == pytest.approx(6.0, abs=1e-10)
    assert lengths[0] == pytest.approx(np.sqrt((4.0**2+4.002**2)/2))
    np.testing.assert_allclose(result.scaled_positions, [[0, 0, 0]])


def test_idealization_preserves_existing_metric_and_origin():
    dataset = SimpleNamespace(
        rotations=[np.eye(3)], translations=[np.zeros(3)]
    )
    cell = np.array([[2.0, 0.3, 0], [0.2, 3.0, 0], [0.1, 0.2, 4.0]])
    positions = np.array([[0.37, 0.28, 0.19]])
    transformed, changed = idealize_cell(cell, positions, [1], dataset)
    np.testing.assert_array_equal(transformed.cell, cell)
    np.testing.assert_allclose(transformed.scaled_positions, positions)
    assert not changed


def test_operation_permutations_shared_with_spacegroup_info():
    from ..spglib import symmetry_site_permutations

    atoms = SPRKKRAtoms(bulk("Fe", "bcc", a=2.87, cubic=True), symmetry=True)
    info = atoms.spacegroup_info
    table = info.kto_kyda_table()
    tags = site_equivalence_tags(atoms)
    expected = symmetry_site_permutations(
        atoms.cell, atoms.get_scaled_positions(),
        tags, info.dataset.rotations, info.dataset.translations,
    )
    np.testing.assert_array_equal(table, expected)



def test_primitive_site_groups_distinguish_same_species():
    # One species tag, two different primitive sites, each repeated by
    # translation. Output order is deliberately reversed.
    cell = np.diag([2.0, 2.0, 2.0])
    positions = np.array([
        [0., 0., 0.], [.5, 0., 0.],
        [.25, .25, .25], [.75, .25, .25],
    ])
    dataset = SimpleNamespace(
        mapping_to_primitive=np.array([0, 0, 1, 1]),
        primitive_lattice=np.diag([1., 2., 2.]),
        transformation_matrix=np.eye(3),
        origin_shift=np.zeros(3),
    )
    transformed = TransformedCell(
        np.diag([1., 2., 2.]),
        np.array([[.5, .25, .25], [0., 0., 0.]]),
        np.array([1, 1]),
    )
    assert primitive_site_groups(
        cell, positions, [1, 1, 1, 1], dataset, transformed,
    ) == [[2, 3], [0, 1]]


def test_primitive_mapping_retains_input_origin_without_idealization():
    cell = np.diag([2., 2., 2.])
    positions = np.array([[0., 0., 0.], [.25, .25, .25]])
    dataset = SimpleNamespace(
        mapping_to_primitive=np.array([0, 1]),
        primitive_lattice=cell,
        transformation_matrix=np.eye(3),
        origin_shift=np.array([.25, 0., 0.]),
    )
    transformed = TransformedCell(
        cell, positions[::-1].copy(),
        np.array([1, 1]),
    )
    assert primitive_site_groups(
        cell, positions, [1, 1], dataset, transformed,
    ) == [[1], [0]]


def test_explicit_kinds_override_automatic_species_partition():
    atoms = bulk("NaCl", "rocksalt", a=5.64)
    normal = site_equivalence_tags(atoms)
    explicit = site_equivalence_tags(atoms, atomic_kinds=[4, 4])
    assert normal[0] != normal[1]
    assert explicit[0] == explicit[1]
    assert symmetry_dataset(
        atoms.cell, atoms.get_scaled_positions(), explicit
    ).equivalent_atoms[0] == 0


@pytest.mark.parametrize("equivalence", ["site-types", "occupation"])
def test_ase_kinds_resolve_composition_without_splitting_equal_sites(equivalence):
    atoms = bulk("Fe", "bcc", a=2.87, cubic=True)
    atoms.set_array("spacegroup_kinds", np.array([0, 1]))
    atoms.info["occupancy"] = {"0": {"Fe": 1.0}, "1": {"Fe": 1.0}}
    tags = site_equivalence_tags(atoms, equivalence=equivalence)
    assert tags[0] == tags[1]
    assert symmetry_dataset(atoms.cell, atoms.get_scaled_positions(), tags).number == 229

    atoms.info["occupancy"]["1"] = {"Fe": 0.5, "Ni": 0.5}
    tags = site_equivalence_tags(atoms, equivalence=equivalence)
    assert tags[0] != tags[1]


def test_previously_assigned_region_sites_refine_explicit_kinds():
    atoms = SPRKKRAtoms(
        bulk("Fe", "bcc", a=2.87, cubic=True), symmetry=False,
    )
    sites = atoms.sites
    tags = site_equivalence_tags(
        atoms, atomic_kinds=[7, 7], assigned_sites=sites,
    )
    assert tags[0] != tags[1]


def test_primitive_correspondence_with_nontrivial_basis_change():
    # Swap the first two primitive basis vectors. No-idealization output
    # retains the Cartesian origin despite the dataset's origin_shift.
    cell = np.diag([2., 3., 4.])
    positions = np.array([[.13, .21, .27], [.43, .51, .17]])
    P = np.array([[0., 1., 0.], [1., 0., 0.], [0., 0., 1.]])
    origin = np.array([.25, .125, 0.])
    standard_cell = np.linalg.solve(P.T, cell)
    transformed = TransformedCell(
        standard_cell,
        (positions[::-1] @ P.T) % 1.,
        np.array([1, 1]),
    )
    dataset = SimpleNamespace(
        mapping_to_primitive=np.array([0, 1]),
        primitive_lattice=cell,
        transformation_matrix=P,
        origin_shift=origin,
    )
    assert primitive_site_groups(
        cell, positions, [1, 1], dataset, transformed
    ) == [[1], [0]]


def test_idealization_symmetrizes_inversion_related_positions():
    cell = np.eye(3) * 4.
    positions = np.array([[.201, .3, .3], [.800, .7, .7]])
    data = SimpleNamespace(
        rotations=np.array([np.eye(3), -np.eye(3)]),
        translations=np.zeros((2, 3)),
    )
    result, changed = idealize_cell(
        cell, positions, [1, 1], data, symprec=.01
    )
    assert changed
    np.testing.assert_allclose(
        (result.scaled_positions[0] + result.scaled_positions[1]) % 1.,
        [0., 0., 0.], atol=1e-12,
    )
    np.testing.assert_allclose(result.cell, cell)


def test_strict_initial_magmoms_are_part_of_site_classification():
    atoms = bulk("Fe", "bcc", a=2.87, cubic=True)
    atoms.set_initial_magnetic_moments([1.0, -1.0])
    strict = site_equivalence_tags(atoms)
    occupation = site_equivalence_tags(atoms, equivalence="occupation")
    assert strict[0] != strict[1]
    assert occupation[0] == occupation[1]


def test_assigned_sites_are_partitioned_by_site_type_not_site_identity():
    atoms = SPRKKRAtoms(bulk("Fe", "bcc", a=2.87, cubic=True))
    sites = atoms.sites
    assert sites[0] is not sites[1]
    assert sites[0].site_type is sites[1].site_type
    tags = site_equivalence_tags(atoms, assigned_sites=sites)
    assert tags[0] == tags[1]


def test_region_memberships_remain_distinct_even_with_shared_site_type():
    atoms = SPRKKRAtoms(bulk("Fe", "bcc", a=2.87, cubic=True))
    sites = atoms.sites
    tags = site_equivalence_tags(
        atoms, assigned_sites=sites,
        region_memberships=[("regionA",), ("regionB",)],
    )
    assert tags[0] != tags[1]


def test_legacy_symmetry_scenarios_use_new_single_dataset_api():
    from ...common.unique_values import UniqueValuesMapping

    def classes(dataset, expected):
        assert UniqueValuesMapping(dataset.equivalent_atoms).is_equivalent_to(expected)

    nacl = bulk("NaCl", "rocksalt", a=5.64)
    classes(
        symmetry_dataset(nacl.cell, nacl.get_scaled_positions(),
                         site_equivalence_tags(nacl)),
        [0, 1],
    )
    classes(
        symmetry_dataset(
            nacl.cell, nacl.get_scaled_positions(),
            site_equivalence_tags(nacl, atomic_kinds=[4, 4]),
        ),
        [0, 0],
    )

    # Two crystallographic positions of the same element.
    nana = bulk("NaNa", "rocksalt", a=5.64)
    classes(
        symmetry_dataset(nana.cell, nana.get_scaled_positions(),
                         site_equivalence_tags(nana)),
        [0, 0],
    )
    classes(
        symmetry_dataset(
            nana.cell, nana.get_scaled_positions(),
            site_equivalence_tags(nana, atomic_kinds=[4, 7]),
        ),
        [0, 1],
    )

    atoms = SPRKKRAtoms(nacl, symmetry=False)
    original = atoms.copy()
    atoms.sites[1].occupation = {"Na": 1.0}
    strict = site_equivalence_tags(atoms)
    relaxed = site_equivalence_tags(atoms, equivalence="occupation")
    assert strict[0] != strict[1]
    assert relaxed[0] == relaxed[1]
    np.testing.assert_allclose(original.positions, atoms.positions)



def test_legacy_dict_dataset_supports_attribute_access(monkeypatch):
    """spglib before 2.5 returns a dict; no minimum spglib version required."""
    import spglib

    sample = {
        "number": 229,
        "equivalent_atoms": np.array([0, 0]),
        "mapping_to_primitive": np.array([0, 0]),
        "primitive_lattice": np.eye(3),
        "transformation_matrix": np.eye(3),
        "origin_shift": np.zeros(3),
        "rotations": np.eye(3, dtype=int)[None, :, :],
        "translations": np.zeros((1, 3)),
    }
    monkeypatch.setattr(spglib, "get_symmetry_dataset", lambda *args, **kw: sample)
    result = symmetry_dataset(np.eye(3), np.zeros((2, 3)), [1, 1])
    assert result.number == 229
    np.testing.assert_array_equal(result.mapping_to_primitive, [0, 0])
    with pytest.raises(AttributeError):
        _ = result.non_existent_value


def test_modern_object_dataset_is_not_wrapped(monkeypatch):
    """Preserve spglib's native dataset when it already has attributes."""
    import spglib

    result = SimpleNamespace(number=225, equivalent_atoms=np.array([0, 1]))
    monkeypatch.setattr(spglib, "get_symmetry_dataset", lambda *args, **kw: result)
    assert symmetry_dataset(np.eye(3), np.zeros((2, 3)), [1, 2]) is result


@pytest.mark.parametrize("api", ["get_symmetry_dataset", "standardize_cell"])
def test_old_spglib_none_failure_becomes_value_error(monkeypatch, api):
    import spglib

    monkeypatch.setattr(spglib, api, lambda *args, **kw: None)
    if api == "get_symmetry_dataset":
        call = lambda: symmetry_dataset(np.eye(3), [[0, 0, 0]], [1])
    else:
        call = lambda: transform_cell(np.eye(3), [[0, 0, 0]], [1])
    with pytest.raises(ValueError, match="spglib"):
        call()


@pytest.mark.parametrize("api", ["get_symmetry_dataset", "standardize_cell"])
def test_new_spglib_exception_becomes_value_error(monkeypatch, api):
    import spglib
    from ..spglib import SpglibError

    def fail(*args, **kwargs):
        raise SpglibError("test failure")

    monkeypatch.setattr(spglib, api, fail)
    if api == "get_symmetry_dataset":
        call = lambda: symmetry_dataset(np.eye(3), [[0, 0, 0]], [1])
    else:
        call = lambda: transform_cell(np.eye(3), [[0, 0, 0]], [1])
    with pytest.raises(ValueError, match="test failure"):
        call()


def test_standardize_preserves_spglib_arrays_without_recasting(monkeypatch):
    """The spglib boundary passes its already-typed output arrays through."""
    import spglib

    lattice = np.eye(3, dtype=np.float64)
    positions = np.array([[0., 0., 0.]], dtype=np.float64)
    kinds = np.array([4], dtype=np.int32)
    monkeypatch.setattr(
        spglib, "standardize_cell",
        lambda *args, **kwargs: (lattice, positions, kinds),
    )

    result = transform_cell(lattice, positions, kinds)
    assert result.cell is lattice
    assert result.scaled_positions is positions
    assert result.tags is kinds


def test_operation_permutations_accept_plain_python_input_lists():
    from ..spglib import symmetry_site_permutations

    table = symmetry_site_permutations(
        cell=[[2., 0., 0.], [0., 2., 0.], [0., 0., 2.]],
        scaled_positions=[[0., 0., 0.], [.5, .5, .5]],
        tags=[1, 1],
        rotations=[np.eye(3, dtype=np.intc)],
        translations=[np.zeros(3)],
    )
    np.testing.assert_array_equal(table, [[0, 1]])
