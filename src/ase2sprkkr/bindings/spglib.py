"""Small, single-purpose spglib binding used by SPR-KKR symmetry workflows."""

from dataclasses import dataclass

import numpy as np
import spglib
from ase.geometry import find_mic

try:
    from spglib.error import SpglibError
except ImportError:
    # Before spglib 2.7 the public API returned None on symmetry-search
    # failure rather than raising SpglibError.
    class SpglibError(Exception):
        pass


class _SpglibCompatibilityWrapper:
    """Provide attribute access to datasets returned as dicts by old spglib."""

    def __init__(self, dataset):
        self.dataset = dataset

    def __getattr__(self, name):
        # Avoid recursion if an unpickled instance is not initialized yet.
        try:
            dataset = self.__dict__["dataset"]
        except KeyError:
            raise AttributeError(name) from None
        try:
            return dataset[name]
        except KeyError:
            raise AttributeError(name) from None

    @classmethod
    def create_if_needed(cls, dataset):
        if dataset is None or hasattr(dataset, "equivalent_atoms"):
            return dataset
        return cls(dataset)


@dataclass(frozen=True)
class TransformedCell:
    """Geometry and species labels in the order supplied by spglib."""

    cell: np.ndarray
    scaled_positions: np.ndarray
    tags: np.ndarray


def symmetry_dataset(cell, scaled_positions, site_tags, *, symprec=1e-5, angle_tolerance=-1.0):
    """Find a space-group dataset for explicitly classified atomic sites."""
    try:
        dataset = spglib.get_symmetry_dataset(
            (cell, scaled_positions, site_tags),
            symprec=symprec, angle_tolerance=angle_tolerance,
        )
    except SpglibError as exc:
        raise ValueError(f"spglib symmetry search failed: {exc}") from exc
    if dataset is None:
        raise ValueError("spglib could not determine the space-group symmetry")
    return _SpglibCompatibilityWrapper.create_if_needed(dataset)


def transform_cell(
    cell, scaled_positions, tags, *, to_primitive=False,
    symprec=1e-5, angle_tolerance=-1.0,
):
    """Standardize without idealization, retaining the spglib species labels.

    Every output tag identifies the source SiteType/occupation class. Site
    metadata transfer is independent of origin, orientation, and atom order.
    """
    try:
        result = spglib.standardize_cell(
            (cell, scaled_positions, tags),
            to_primitive=to_primitive, no_idealize=True,
            symprec=symprec, angle_tolerance=angle_tolerance,
        )
    except SpglibError as exc:
        raise ValueError(f"spglib cell standardization failed: {exc}") from exc
    if result is None:
        raise ValueError("spglib could not standardize the unit cell")
    # spglib already returns NumPy arrays with the required numeric dtypes.
    return TransformedCell(*result)


def primitive_site_groups(cell, scaled_positions, tags, dataset, transformed, *, symprec=1e-5):
    """Relate each output atom to its input primitive translation class.

    mapping_to_primitive distinguishes inequivalent primitive sites, unlike
    the chemical-species tags. We align input and standardized output frames
    using spglib's P and origin shift, without assuming that numeric primitive
    IDs of independently generated cells have identical meaning.
    """
    # Only normalize caller-provided inputs; dataset fields are ndarrays.
    cell = np.asarray(cell, dtype=float)
    positions = np.asarray(scaled_positions, dtype=float)
    tags = np.asarray(tags)
    primitive_ids = dataset.mapping_to_primitive
    if len(primitive_ids) != len(tags):
        raise ValueError("Invalid spglib mapping_to_primitive length")
    classes = {}
    for i, primitive_id in enumerate(primitive_ids):
        classes.setdefault(int(primitive_id), []).append(i)

    matrix = dataset.transformation_matrix
    standard_cell = np.linalg.solve(matrix.T, cell)
    cartesian_shift = dataset.origin_shift @ standard_cell
    input_cart = positions @ cell + cartesian_shift
    output_cart = transformed.scaled_positions @ transformed.cell
    primitive_cell = dataset.primitive_lattice

    by_tag = {}
    for group in classes.values():
        values = {int(tags[i]) for i in group}
        if len(values) != 1:
            raise ValueError("A primitive site contains multiple species")
        by_tag.setdefault(values.pop(), []).append(group)

    output_groups = []
    tolerance = max(5 * symprec, 1e-7)
    for output_index, (point, tag) in enumerate(zip(output_cart, transformed.tags)):
        possible = by_tag.get(int(tag), [])
        if not possible:
            raise ValueError(f"Unknown output species tag {tag}")
        representatives = [group[0] for group in possible]
        _, distances = find_mic(
            input_cart[representatives] - point, cell=primitive_cell, pbc=True
        )
        best = int(np.argmin(distances))
        if distances[best] > tolerance:
            raise ValueError(
                f"Cannot map transformed atom {output_index} to a primitive site: "
                f"nearest distance {distances[best]:g} Å"
            )
        if np.count_nonzero(distances <= tolerance) != 1:
            raise ValueError(f"Ambiguous primitive class for transformed atom {output_index}")
        output_groups.append(possible[best])
    return output_groups


def symmetry_site_permutations(
    cell, scaled_positions, tags, rotations, translations, *, symprec=1e-5
):
    """Return the atomic permutation for each space-group operation.

    The same periodic, species-aware matching is used by both metric/site
    idealization and SpacegroupInfo.kto_kyda_table.
    """
    cell = np.asarray(cell, dtype=float)
    positions = np.asarray(scaled_positions, dtype=float)
    tags = np.asarray(tags)
    by_tag = {tag: np.flatnonzero(tags == tag) for tag in np.unique(tags)}
    table = np.empty((len(rotations), len(positions)), dtype=int)
    tolerance = max(2 * symprec, 1e-8)
    for g, (rotation, translation) in enumerate(zip(rotations, translations)):
        for i, (position, tag) in enumerate(zip(positions, tags)):
            transformed = position @ rotation.T + translation
            candidates = by_tag[tag]
            _, distances = find_mic(
                (positions[candidates] - transformed) @ cell, cell=cell, pbc=True
            )
            j = int(np.argmin(distances))
            if distances[j] > tolerance:
                raise ValueError(
                    f"Cannot map operation {g}, site {i} within {tolerance:g} Å"
                )
            table[g, i] = candidates[j]
        if len(set(table[g])) != len(positions):
            raise ValueError(f"Operation {g} does not permute all sites")
    return table


def idealize_cell(cell, scaled_positions, tags, dataset, *, symprec=1e-5):
    """Symmetrize lattice metric and atomic positions in the ORIGINAL basis.

    The metric is projected onto the point-group-invariant subspace.
    For each atomic image, the inverse symmetry operation maps it back
    near its source site; minimum-image displacements are then averaged.
    This keeps the incoming basis setting, atom count, and ordering.
    """
    cell = np.asarray(cell, dtype=float)
    positions = np.asarray(scaled_positions, dtype=float)
    tags = np.asarray(tags)
    rotations = dataset.rotations
    translations = dataset.translations

    metric = cell @ cell.T
    sym_metric = np.mean(
        [rotation.T @ metric @ rotation for rotation in rotations], axis=0
    )
    sym_metric = 0.5 * (sym_metric + sym_metric.T)

    def symmetric_root(matrix, inverse=False):
        eigenvalues, eigenvectors = np.linalg.eigh(matrix)
        if np.any(eigenvalues <= 0):
            raise ValueError("A positive-definite lattice metric is required")
        powers = eigenvalues ** (-0.5 if inverse else 0.5)
        return (eigenvectors * powers) @ eigenvectors.T

    new_cell = symmetric_root(sym_metric) @ symmetric_root(metric, inverse=True) @ cell
    inverse_cell = np.linalg.inv(cell)

    delta_sum = np.zeros_like(positions)
    permutations = symmetry_site_permutations(
        cell, positions, tags, rotations, translations, symprec=symprec
    )
    for rotation, translation, mapping in zip(rotations, translations, permutations):
        inverse_rotation = np.linalg.inv(rotation.T)
        for i, j in enumerate(mapping):
            back = (positions[j] - translation) @ inverse_rotation
            correction, _ = find_mic(
                (back - positions[i]) @ cell, cell=cell, pbc=True
            )
            delta_sum[i] += correction @ inverse_cell

    ideal_positions = (positions + delta_sum / len(rotations)) % 1.0
    metric_change = not np.allclose(
        metric, sym_metric, rtol=1e-10, atol=1e-10
    )
    # Compare displacement, not modulo-wrapped coordinates.
    displacement = np.linalg.norm((delta_sum / len(rotations)) @ cell, axis=1)
    geometry_changed = metric_change or bool(np.any(displacement > 1e-8))
    return TransformedCell(new_cell, ideal_positions, tags.copy()), geometry_changed
