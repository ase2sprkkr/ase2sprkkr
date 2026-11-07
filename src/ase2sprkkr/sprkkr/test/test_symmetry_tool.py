"""Crystal symmetry, site ownership, ASE-array policy, and CLI tests."""

import argparse
import warnings
from unittest.mock import patch

import numpy as np
import pytest
from ase.build import bulk
from ase.io import write, read

from ..sprkkr_atoms import SPRKKRAtoms
from ..symmetry import SymmetryError, analyze_symmetry
from ...tools.commands import symmetry as symmetry_command


def test_primitive_bcc_from_plain_ase():
    original = bulk("Fe", "bcc", a=2.87, cubic=True)
    before = original.copy()
    analysis = analyze_symmetry(original)
    assert analysis.dataset.number == 229
    primitive = analysis.find_primitive(original)
    assert isinstance(primitive, SPRKKRAtoms)
    assert len(primitive) == 1
    assert len(original) == 2
    np.testing.assert_allclose(original.positions, before.positions)
    assert primitive.sites[0].site_type.atoms is primitive


def test_nacl_primitive():
    original = bulk("NaCl", "rocksalt", a=5.64, cubic=True)
    analysis = analyze_symmetry(original)
    primitive = analysis.find_primitive(original)
    assert analysis.dataset.number == 225
    assert len(original) == 8
    assert len(primitive) == 2
    assert set(primitive.get_chemical_symbols()) == {"Na", "Cl"}


def test_standardize_fcc_from_primitive():
    original = bulk("Cu", "fcc", a=3.6)
    conventional = analyze_symmetry(original).standardize(original)
    assert len(conventional) == 4
    assert conventional.sites[0].site_type.atoms is conventional


def test_strict_site_types_vs_occupation():
    original = SPRKKRAtoms(
        bulk("Fe", "bcc", a=2.87).repeat((2, 1, 1)), symmetry=False,
    )
    original.sites
    assert original.sites[0].site_type is not original.sites[1].site_type
    strict = analyze_symmetry(original).find_primitive(original)
    relaxed = analyze_symmetry(original, equivalence="occupation").find_primitive(original)
    assert len(strict) == 2
    assert len(relaxed) == 1


def test_occupation_copy_independence():
    original = SPRKKRAtoms(
        bulk("Fe", "bcc", a=2.87).repeat((2, 1, 1)), symmetry=False,
    )
    original.sites[0].occupation = {"Fe": 0.5, "Ni": 0.5}
    original.sites[1].occupation = {"Ni": 0.5, "Fe": 0.5}
    old_type = original.sites[0].site_type
    result = analyze_symmetry(original, equivalence="occupation").find_primitive(original)
    assert len(result) == 1
    assert result.sites[0].occupation.to_dict(True) == {"Fe": 0.5, "Ni": 0.5}
    assert result.sites[0].site_type is not old_type
    assert result.sites[0].site_type.atoms is result
    result.sites[0].occupation = {"Fe": 1.0}
    assert original.sites[0].occupation.to_dict(True) == {"Fe": 0.5, "Ni": 0.5}


def test_site_type_sharing_in_expanded_cell():
    atoms = SPRKKRAtoms(bulk("Cu", "fcc", a=3.6), symmetry=False)
    output = analyze_symmetry(atoms).standardize(atoms)
    assert len(output) == 4
    assert all(site.site_type is output.sites[0].site_type for site in output.sites)
    assert output.sites[0].site_type is not atoms.sites[0].site_type


def test_ase_array_policies_are_explicit():
    atoms = bulk("Fe", "bcc", a=2.87, cubic=True)
    atoms.set_array("custom", np.asarray([1.0, 2.0]))
    analysis = analyze_symmetry(atoms)

    with pytest.raises(SymmetryError, match="custom"):
        analysis.find_primitive(atoms)

    with pytest.warns(RuntimeWarning, match="first representative"):
        first = analysis.find_primitive(atoms, ase_arrays="choose-first")
    assert list(first.arrays["custom"]) == [1.0]

    with pytest.warns(RuntimeWarning, match="Discarded ASE"):
        discarded = analysis.find_primitive(atoms, ase_arrays="discard")
    assert "custom" not in discarded.arrays

    atoms.set_array("custom", np.asarray([7.0, 7.0]))
    equal = analyze_symmetry(atoms).find_primitive(atoms)
    assert list(equal.arrays["custom"]) == [7.0]


def test_ase_array_is_replicated_on_cell_expansion():
    atoms = bulk("Cu", "fcc", a=3.6)
    atoms.set_array("custom", np.asarray([12]))
    standard = analyze_symmetry(atoms).standardize(atoms)
    assert len(standard) == 4
    np.testing.assert_array_equal(standard.arrays["custom"], [12, 12, 12, 12])


def test_idealize_alone_preserves_basis_order_and_extra_arrays():
    atoms = bulk("Fe", "bcc", a=2.87, cubic=True)
    atoms.set_array("custom", np.asarray([4, 5]))
    before = atoms.copy()
    output = analyze_symmetry(atoms).idealize(atoms)

    assert len(output) == len(atoms)
    np.testing.assert_array_equal(output.arrays["custom"], [4, 5])
    np.testing.assert_array_equal(atoms.arrays["custom"], before.arrays["custom"])
    np.testing.assert_allclose(atoms.cell, before.cell)


def test_idealize_warns_for_electronic_data_only_when_geometry_changes():
    atoms = SPRKKRAtoms(bulk("Fe", "bcc", a=2.87), symmetry=False)
    atom_type = atoms.sites[0].site_type.atomic_type(0)
    atom_type.potential = np.ones((2, len(atom_type.mesh)))
    analysis = analyze_symmetry(atoms)
    from ...bindings.spglib import TransformedCell
    source = atoms
    fake = TransformedCell(
        source.cell.array, source.get_scaled_positions(wrap=True), analysis.site_tags,
    )

    with patch("ase2sprkkr.sprkkr.symmetry.idealize_cell", return_value=(fake, True)):
        with pytest.warns(RuntimeWarning, match="new SCF"):
            result = analysis.idealize(atoms)
    assert result.sites[0].site_type.atomic_type(0).potential is not None

    with patch("ase2sprkkr.sprkkr.symmetry.idealize_cell", return_value=(fake, False)):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            analysis.idealize(atoms)
    assert not any("new SCF" in str(w.message) for w in caught)


def test_idealize_and_standardization_are_independent():
    atoms = bulk("Cu", "fcc", a=3.6)
    analysis = analyze_symmetry(atoms)
    ideal = analysis.idealize(atoms)
    standard = analysis.standardize(atoms)
    both = analysis.standardize(atoms, idealize=True)
    assert len(ideal) == len(atoms)
    assert len(standard) == len(both) == 4


def test_nonperiodic_is_rejected():
    atoms = bulk("Fe", "bcc", a=2.87)
    atoms.pbc = (True, True, False)
    with pytest.raises(SymmetryError, match="periodic"):
        analyze_symmetry(atoms)


def test_command_roundtrip_and_idealize_only(tmp_path, capsys):
    source = tmp_path / "source.cif"
    primitive = tmp_path / "primitive.cif"
    ideal = tmp_path / "ideal.cif"
    write(source, bulk("Fe", "bcc", a=2.87, cubic=True))
    parser = argparse.ArgumentParser()
    symmetry_command.parser(parser)
    symmetry_command.run(parser.parse_args([str(source), "--primitive", "-o", str(primitive)]))
    symmetry_command.run(parser.parse_args([str(source), "--idealize", "-o", str(ideal)]))
    assert len(read(primitive)) == 1
    assert len(read(ideal)) == 2
    assert "Space group" in capsys.readouterr().out


def test_does_not_initialize_caller_sites():
    atoms = SPRKKRAtoms(bulk("Fe", "bcc", a=2.87), symmetry=False)
    assert not atoms.are_sites_inited()
    info = atoms.spacegroup_info
    analysis = analyze_symmetry(atoms)
    assert analysis.dataset.number == 229
    assert not atoms.are_sites_inited()
    assert atoms.spacegroup_info is info
    assert atoms.spacegroup_info._dataset is None


def test_distinct_primitive_sites_of_same_species_keep_distinct_ase_values():
    # The Fe atoms are not related by a primitive translation, although
    # they have the same spglib species tag.
    atoms = bulk("Fe", "bcc", a=2.87, cubic=True)
    positions = atoms.get_scaled_positions()
    positions[1] += [0.13, 0.07, 0.0]
    atoms.set_scaled_positions(positions)
    atoms.set_array("local", np.asarray([10, 20]))
    result = analyze_symmetry(atoms).find_primitive(atoms)
    assert len(result) == 2
    assert sorted(result.arrays["local"].tolist()) == [10, 20]


def test_ase_equality_only_within_translation_classes():
    from ase import Atoms

    # Two inequivalent primitive sites, each repeated along a.
    atoms = Atoms(
        "Fe4",
        scaled_positions=[
            [0., 0., 0.], [.25, .2, .1],
            [.5, 0., 0.], [.75, .2, .1],
        ],
        cell=[4., 3., 3.], pbc=True,
    )
    atoms.set_array("local", np.array([10, 20, 10, 20]))
    primitive = analyze_symmetry(atoms).find_primitive(atoms)
    assert len(primitive) == 2
    assert sorted(primitive.arrays["local"].tolist()) == [10, 20]
    atoms.arrays["local"][2] = 30
    with pytest.raises(SymmetryError, match="local"):
        analyze_symmetry(atoms).find_primitive(atoms)


def test_ase_array_policy_can_choose_initial_magmoms():
    atoms = bulk("Fe", "bcc", a=2.87, cubic=True)
    atoms.set_initial_magnetic_moments([1.0, 2.0])
    analysis = analyze_symmetry(atoms, equivalence="occupation")
    with pytest.raises(SymmetryError, match="initial_magmoms"):
        analysis.find_primitive(atoms)
    with pytest.warns(RuntimeWarning, match="initial_magmoms"):
        result = analysis.find_primitive(atoms, ase_arrays="choose-first")
    assert len(result) == 1
    assert result.get_initial_magnetic_moments()[0] in (1.0, 2.0)


def test_primitive_conversion_copies_radial_potential_independently():
    atoms = SPRKKRAtoms(
        bulk("Fe", "bcc", a=2.87, cubic=True), symmetry=True,
    )
    source_type = atoms.sites[0].site_type
    atype = source_type.atomic_type(0)
    atype.potential = np.ones((2, len(atype.mesh)))

    primitive = analyze_symmetry(atoms).find_primitive(atoms)
    target_type = primitive.sites[0].site_type
    target_atomic = target_type.atomic_type(0)
    assert target_type is not source_type
    assert target_type.atoms is primitive
    assert target_atomic.potential is not atype.potential
    np.testing.assert_allclose(
        target_atomic.potential.raw_value,
        atype.potential.raw_value,
    )


def test_strict_initial_magmoms_prevent_reduction_without_data_loss():
    atoms = bulk("Fe", "bcc", a=2.87, cubic=True)
    atoms.set_initial_magnetic_moments([2.0, -2.0])
    result = analyze_symmetry(atoms).find_primitive(atoms)
    assert len(result) == 2
    np.testing.assert_allclose(
        np.sort(result.get_initial_magnetic_moments()), [-2.0, 2.0]
    )


@pytest.mark.parametrize("structure, primitive_count, conventional_count", [
    ("NaCl", 2, 8),
    ("C", 2, 8),
])
def test_real_spglib_primitive_and_conventional_with_shifted_origin(
    structure, primitive_count, conventional_count,
):
    if structure == "NaCl":
        atoms = bulk("NaCl", "rocksalt", a=5.64)
    else:
        atoms = bulk("C", "diamond", a=3.57)
    atoms.set_scaled_positions(
        atoms.get_scaled_positions() + [.173, .221, .091]
    )
    analysis = analyze_symmetry(atoms)
    primitive = analysis.find_primitive(atoms)
    conventional = analysis.standardize(atoms)
    assert len(primitive) == primitive_count
    assert len(conventional) == conventional_count
    np.testing.assert_allclose(
        abs(np.linalg.det(primitive.cell)),
        abs(np.linalg.det(atoms.cell)),
        rtol=1e-6,
    )
    np.testing.assert_allclose(
        abs(np.linalg.det(conventional.cell)),
        4 * abs(np.linalg.det(primitive.cell)),
        rtol=1e-6,
    )


def test_real_spglib_two_diamond_primitive_sites_same_species():
    atoms = bulk("C", "diamond", a=3.57)
    atoms.set_array("local", np.array([7, 13]))
    conventional = analyze_symmetry(atoms).standardize(atoms)
    assert len(conventional) == 8
    np.testing.assert_array_equal(
        np.sort(conventional.arrays["local"]), [7] * 4 + [13] * 4
    )


def test_real_spglib_sheared_primitive_and_shifted_origin():
    atoms = bulk("C", "diamond", a=3.57)
    cartesian = atoms.get_positions().copy()
    matrix = np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
    atoms.set_cell(matrix @ atoms.cell.array, scale_atoms=False)
    np.testing.assert_allclose(atoms.positions, cartesian)
    atoms.set_scaled_positions(
        atoms.get_scaled_positions() + [.14, .23, .09]
    )
    atoms.set_array("local", np.array([5, 9]))
    transformed = analyze_symmetry(atoms).standardize(atoms)
    assert len(transformed) == 8
    np.testing.assert_array_equal(
        np.sort(transformed.arrays["local"]), [5] * 4 + [9] * 4
    )


def test_potential_serialization_roundtrip_after_cell_conversion():
    import io
    from ...potentials.potentials import Potential

    atoms = bulk("NaCl", "rocksalt", a=5.64, cubic=True)
    transformed = analyze_symmetry(atoms).find_primitive(atoms)
    buffer = io.StringIO()
    Potential.from_atoms(transformed).save_to_file(buffer)
    buffer.seek(0)
    restored = Potential.from_file(buffer).atoms
    assert len(restored) == 2
    np.testing.assert_allclose(restored.cell, transformed.cell)
    np.testing.assert_allclose(restored.positions, transformed.positions)
    assert sorted(restored.get_chemical_symbols()) == ["Cl", "Na"]


def test_conversion_preserves_nonstructural_ase_info():
    atoms = bulk("Fe", "bcc", a=2.87, cubic=True)
    atoms.info["label"] = {"sample": ["alpha"]}
    result = analyze_symmetry(atoms).find_primitive(atoms)
    assert result.info["label"] == {"sample": ["alpha"]}
    result.info["label"]["sample"].append("beta")
    assert atoms.info["label"] == {"sample": ["alpha"]}


def test_analysis_never_copies_or_promotes_plain_ase_atoms(monkeypatch):
    from ase import Atoms

    atoms = bulk("Fe", "bcc", a=2.87, cubic=True)

    def forbidden(*args, **kwargs):
        raise AssertionError("symmetry analysis must not copy or promote Atoms")

    monkeypatch.setattr(atoms, "copy", forbidden)
    monkeypatch.setattr(SPRKKRAtoms, "promote_ase_atoms", forbidden)
    analysis = analyze_symmetry(atoms)
    assert analysis.dataset.number == 229
    assert all(not isinstance(v, Atoms) for v in vars(analysis).values())
    assert type(atoms) is Atoms


def test_analysis_keeps_uninitialized_sprkkr_atoms_unmodified(monkeypatch):
    atoms = SPRKKRAtoms(bulk("Fe", "bcc", a=2.87), symmetry=False)
    assert not atoms.are_sites_inited()
    with patch.object(atoms, "copy", side_effect=AssertionError("copy was called")):
        analysis = analyze_symmetry(atoms)
    assert analysis.dataset.number == 229
    assert not atoms.are_sites_inited()


def test_analysis_with_partial_sites_never_repairs_source_in_place():
    atoms = SPRKKRAtoms(bulk("Fe", "bcc", a=2.87, cubic=True), symmetry=False)
    original = atoms.sites.copy()
    partial = original.copy()
    partial[0] = 0
    atoms.set_array("sprkkr_sites_data", partial)
    analysis = analyze_symmetry(atoms)
    assert analysis.dataset is not None
    assert atoms.arrays["sprkkr_sites_data"][0] == 0
    assert atoms.arrays["sprkkr_sites_data"][1] is original[1]


def test_analysis_rejects_atoms_changed_after_dataset_was_computed():
    atoms = bulk("Fe", "bcc", a=2.87, cubic=True)
    analysis = analyze_symmetry(atoms)
    atoms.positions[0, 0] += 0.1
    with pytest.raises(SymmetryError, match="Geometry changed"):
        analysis.find_primitive(atoms)


def test_conversion_without_source_site_initialization():
    atoms = SPRKKRAtoms(bulk("Fe", "bcc", a=2.87, cubic=True), symmetry=True)
    assert not atoms.are_sites_inited()
    result = analyze_symmetry(atoms).find_primitive(atoms)
    assert len(result) == 1
    assert not atoms.are_sites_inited()
    assert result.sites[0].site_type.atoms is result


def test_plain_ase_idealization_keeps_plain_ase_type():
    from ase import Atoms

    atoms = bulk("Fe", "bcc", a=2.87, cubic=True)
    result = analyze_symmetry(atoms).idealize(atoms)
    assert type(result) is Atoms
    assert result is not atoms


def test_site_classification_change_invalidates_cached_analysis():
    atoms = bulk("Fe", "bcc", a=2.87, cubic=True)
    analysis = analyze_symmetry(atoms)
    atoms.set_initial_magnetic_moments([1.0, -1.0])
    with pytest.raises(SymmetryError, match="Site classification changed"):
        analysis.find_primitive(atoms)
