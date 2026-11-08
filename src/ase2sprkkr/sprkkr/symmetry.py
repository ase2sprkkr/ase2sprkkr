"""Crystal symmetry analysis without retaining or modifying input Atoms."""

from __future__ import annotations

import copy
import warnings

import numpy as np
from ase import Atoms

from ..bindings.spglib import (
    idealize_cell, primitive_site_groups, symmetry_dataset, transform_cell,
)
from .occupations import Occupation
from .site_equivalence import site_equivalence_tags, refine_partial_site_tags
from .sites import Site, SiteType
from .sprkkr_atoms import SPRKKRAtoms


class SymmetryError(ValueError):
    """Unable to analyse or transform a crystal."""


def _validate(atoms):
    if not isinstance(atoms, Atoms):
        raise TypeError("Expected an ASE Atoms object")
    if len(atoms) == 0 or not np.all(atoms.pbc):
        raise SymmetryError("Expected a nonempty, fully periodic 3D structure")
    cell = atoms.cell.array
    if not np.isfinite(cell).all() or np.linalg.matrix_rank(cell) != 3:
        raise SymmetryError("The cell must have three independent finite vectors")
    if isinstance(atoms, SPRKKRAtoms) and atoms.regions:
        raise SymmetryError("Structures with AtomsRegion are not supported")


def _classify_sites(atoms, equivalence, cell, positions, *, symprec, angle_tolerance):
    """Determine spglib species tags without creating or modifying Site objects."""
    tags = site_equivalence_tags(atoms, equivalence=equivalence)
    source_sites = atoms.arrays.get(SPRKKRAtoms.sites_array_name)
    if (
        equivalence == "site-types" and source_sites is not None
        and any(not site for site in source_sites)
    ):
        preliminary = symmetry_dataset(
            cell, positions, tags,
            symprec=symprec, angle_tolerance=angle_tolerance,
        )
        tags = refine_partial_site_tags(tags, preliminary, source_sites)
    return tags


def _electronic_data_present(atoms):
    sites = atoms.arrays.get(SPRKKRAtoms.sites_array_name)
    if sites is None:
        return False
    return any(
        atomic_type.potential is not None
        or atomic_type.charge is not None
        or atomic_type.moments is not None
        for site in sites if site
        for atomic_type in site.site_type.occupation.atomic_types()
    )


def _warn_electronic_data():
    warnings.warn(
        "Idealization changed the geometry; copied SPR-KKR electronic data "
        "may require a new SCF calculation.",
        RuntimeWarning, stacklevel=3,
    )


def _make_sites(atoms, result, source_indices, dataset):
    """Copy existing SPR-KKR sites, or create new sites for ordinary ASE atoms.

    This function only reads the source. It never initializes source.sites.
    Symmetry-equivalent output atoms share their new SiteType.
    """
    source_sites = atoms.arrays.get(SPRKKRAtoms.sites_array_name)
    known_by_orbit = {}
    if source_sites is not None:
        for i, site in enumerate(source_sites):
            if site:
                known_by_orbit.setdefault(int(dataset.equivalent_atoms[i]), site)

    site_types = {}
    sites = []
    for index in source_indices:
        old_site = (
            source_sites[index]
            if source_sites is not None and source_sites[index]
            else known_by_orbit.get(int(dataset.equivalent_atoms[index]))
        )
        if old_site:
            key = ("site-type", old_site.site_type)
            if key not in site_types:
                site_types[key] = old_site.site_type.copy(atoms=result)
            sites.append(old_site.copy(site_type=site_types[key]))
        else:
            key = ("symmetry-orbit", int(dataset.equivalent_atoms[index]))
            if key not in site_types:
                occupation = Occupation.for_atom_from_ase_atoms(atoms, int(index))
                site_types[key] = SiteType(result, occupation)
            sites.append(Site(site_types[key]))
    return np.asarray(sites, dtype=object)


class SymmetryAnalysis:
    """Symmetry dataset plus small geometry arrays, never a retained Atoms.

    Calling methods on atoms with changed positions, lattice, or site classes
    raises SymmetryError; recompute the analysis for the new structure.
    """

    def __init__(
        self, cell, positions, site_tags, dataset, *,
        equivalence, symprec, angle_tolerance,
    ):
        self._cell = cell.copy()
        self._positions = positions.copy()
        self.site_tags = site_tags.copy()
        self.dataset = dataset
        self.equivalence = equivalence
        self.symprec = symprec
        self.angle_tolerance = angle_tolerance

    def _check_atoms(self, atoms):
        _validate(atoms)
        if (
            len(atoms) != len(self.site_tags)
            or not np.array_equal(atoms.cell.array, self._cell)
            or not np.allclose(
                atoms.get_scaled_positions(wrap=True), self._positions,
                rtol=0, atol=1e-12,
            )
        ):
            raise SymmetryError(
                "Geometry changed since analysis; call analyze_symmetry(atoms) again"
            )
        try:
            tags = _classify_sites(
                atoms, self.equivalence, self._cell, self._positions,
                symprec=self.symprec, angle_tolerance=self.angle_tolerance,
            )
        except ValueError as exc:
            raise SymmetryError(str(exc)) from exc
        if not np.array_equal(tags, self.site_tags):
            raise SymmetryError(
                "Site classification changed since analysis; "
                "call analyze_symmetry(atoms) again"
            )

    def idealize(self, atoms):
        """Idealize a copy in the original basis without modifying atoms."""
        self._check_atoms(atoms)
        if atoms.constraints:
            raise SymmetryError("Idealization with ASE constraints is not supported")
        geometry, changed = idealize_cell(
            self._cell, self._positions, self.site_tags,
            self.dataset, symprec=self.symprec,
        )
        result = atoms.copy()
        result.set_cell(geometry.cell, scale_atoms=False)
        result.set_scaled_positions(geometry.scaled_positions)
        if isinstance(result, SPRKKRAtoms):
            result.spacegroup_info._dataset = None
        if changed and _electronic_data_present(atoms):
            _warn_electronic_data()
        return result

    def find_primitive(self, atoms, *, idealize=False, ase_arrays="require-equal"):
        """Build a primitive cell from explicitly supplied atoms."""
        return self._converted(
            atoms, to_primitive=True, idealize=idealize, ase_arrays=ase_arrays,
        )

    def standardize(self, atoms, *, idealize=False, ase_arrays="require-equal"):
        """Build a standard conventional cell from explicitly supplied atoms."""
        return self._converted(
            atoms, to_primitive=False, idealize=idealize, ase_arrays=ase_arrays,
        )

    def _converted(self, atoms, *, to_primitive, idealize, ase_arrays):
        self._check_atoms(atoms)
        if idealize:
            idealized = self.idealize(atoms)
            updated = analyze_symmetry(
                idealized, equivalence=self.equivalence,
                symprec=self.symprec, angle_tolerance=self.angle_tolerance,
            )
            return updated._converted(
                idealized, to_primitive=to_primitive,
                idealize=False, ase_arrays=ase_arrays,
            )

        if ase_arrays not in ("require-equal", "discard", "choose-first"):
            raise ValueError("ase_arrays must be require-equal, discard, or choose-first")
        if atoms.constraints:
            raise SymmetryError("Cell conversion with ASE constraints is not supported")
        try:
            transformed = transform_cell(
                self._cell, self._positions, self.site_tags,
                to_primitive=to_primitive,
                symprec=self.symprec, angle_tolerance=self.angle_tolerance,
            )
            groups = primitive_site_groups(
                self._cell, self._positions, self.site_tags,
                self.dataset, transformed, symprec=self.symprec,
            )
        except ValueError as exc:
            raise SymmetryError(str(exc)) from exc

        source_sites = atoms.arrays.get(SPRKKRAtoms.sites_array_name)
        source_indices = np.array([
            next(
                (i for i in group if source_sites is not None and source_sites[i]),
                group[0],
            )
            for group in groups
        ], dtype=int)
        result = SPRKKRAtoms(
            numbers=atoms.numbers[source_indices],
            scaled_positions=transformed.scaled_positions,
            cell=transformed.cell, pbc=True,
            symmetry=atoms.symmetry if isinstance(atoms, SPRKKRAtoms) else True,
        )
        result.set_sites(
            _make_sites(atoms, result, source_indices, self.dataset),
            spacegroup_info=True,
        )
        result.spacegroup_info._dataset = None
        result.info.update({
            key: copy.deepcopy(value)
            for key, value in atoms.info.items()
            if key not in ("spacegroup_info", "occupancy")
        })

        reserved = {"numbers", "positions", "sprkkr_sites_data", "spacegroup_kinds"}
        names = set(atoms.arrays) - reserved
        if ase_arrays == "discard":
            if names:
                warnings.warn(
                    "Discarded ASE arrays: " + ", ".join(sorted(names)),
                    RuntimeWarning, stacklevel=2,
                )
        else:
            for name in sorted(names):
                array = atoms.arrays[name]
                if ase_arrays == "require-equal":
                    for indices in groups:
                        if not all(
                            np.array_equal(array[i], array[indices[0]])
                            for i in indices[1:]
                        ):
                            raise SymmetryError(
                                f"ASE array {name!r} differs within a primitive "
                                "translation class; use ase_arrays='choose-first' "
                                "or 'discard'."
                            )
                result.set_array(name, array[source_indices].copy())

            if ase_arrays == "choose-first":
                lossy = [
                    name for name in names
                    if any(
                        not all(
                            np.array_equal(
                                atoms.arrays[name][i],
                                atoms.arrays[name][indices[0]],
                            )
                            for i in indices[1:]
                        )
                        for indices in groups
                    )
                ]
                if lossy:
                    warnings.warn(
                        "Chose the first representative for differing ASE arrays: "
                        + ", ".join(sorted(lossy)),
                        RuntimeWarning, stacklevel=2,
                    )

        if self.equivalence == "occupation" and source_sites is not None:
            original = {site.site_type for site in source_sites if site}
            retained = {
                source_sites[int(i)].site_type
                for i in source_indices if source_sites[int(i)]
            }
            if original - retained:
                warnings.warn(
                    f"Occupation equivalence discarded {len(original - retained)} "
                    "source SiteType(s) and their SPR-KKR data.",
                    RuntimeWarning, stacklevel=2,
                )
        return result


def analyze_symmetry(
    atoms: Atoms, *, equivalence="site-types",
    symprec=1e-5, angle_tolerance=-1.0,
) -> SymmetryAnalysis:
    """Analyse the supplied atoms read-only, without copying or retaining it."""
    _validate(atoms)
    if not np.isfinite(symprec) or symprec <= 0:
        raise ValueError("symprec must be a positive finite length")
    cell = atoms.cell.array
    positions = atoms.get_scaled_positions(wrap=True)
    try:
        tags = _classify_sites(
            atoms, equivalence, cell, positions,
            symprec=symprec, angle_tolerance=angle_tolerance,
        )
        dataset = symmetry_dataset(
            cell, positions, tags,
            symprec=symprec, angle_tolerance=angle_tolerance,
        )
    except ValueError as exc:
        raise SymmetryError(str(exc)) from exc
    return SymmetryAnalysis(
        cell, positions, tags, dataset, equivalence=equivalence,
        symprec=symprec, angle_tolerance=angle_tolerance,
    )
