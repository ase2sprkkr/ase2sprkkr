"""Definition of the SPR-KKR shape-function (``.sfn``) file."""

from collections.abc import Sequence
import math

import numpy as np

from ..output_files import OutputFile
from ..output_files_definitions import (
    OutputFileDefinition,
    OutputFileSectionDefinition,
    OutputFileValueDefinition,
)
from ...common.configuration_definitions import gather
from ...common.dependencies import DependentValue
from ...common.generated_configuration_definitions import GeneratedValueDefinition
from ...common.grammar import line_end
from ...common.grammar_types import Integer, NumpyArray, Real
from ...gui.plot import Multiplot, single_plot


class _SFNValueDefinition(OutputFileValueDefinition):

    """An output value without the leading blank used by standard outputs."""
    prefix = ""


class _SFNSectionDefinition(OutputFileSectionDefinition):
    """A section whose surrounding grammar owns its trailing newline."""

    write_last_delimiter = False


class _ArrayItem:
    """Select an integer item from an array dependency."""

    def __init__(self, index):
        self.index = index

    def __call__(self, value):
        return int(value[self.index])

    @classmethod
    def as_dependent_value(cls, source, *args):
        return DependentValue(source, cls(*args))


class _ArrayShape(_ArrayItem):
    """Build an array shape from one item and fixed trailing dimensions."""

    def __init__(self, index, *trailing):
        super().__init__(index)
        self.trailing = trailing

    def __call__(self, value):
        return (super().__call__(value),) + self.trailing


def _fortran_record(*items):
    """Gather scalar options into one fixed-width record without separators."""

    items = gather(*items)
    items[0].output_definition.value_delimiter = ""
    return items


class ShapeFunctionMesh:
    """Generated combined view of one radial and one Voronoi mesh block."""

    def __init__(self, output, index):
        self._output = output
        self._index = index

    @property
    def radial(self):
        return self._output.RADIAL_MESHES[self._index]

    @property
    def voronoi(self):
        return self._output.VORONOI_MESHES[self._index]

    @property
    def idx(self):
        return self.voronoi.IM()

    @property
    def npan(self):
        return self.radial.NPAN()

    @property
    def nr(self):
        return self.radial.NR()

    @property
    def jrcut(self):
        return self.radial.JRCUT()

    @property
    def radial_mesh(self):
        """The native, unscaled ``(r, dr)`` array stored in the file."""

        return self.radial.RMESH()

    @property
    def rmesh(self):
        return self.radial_mesh[:, 0] * self._output.alat

    @property
    def drmesh(self):
        return self.radial_mesh[:, 1] * self._output.alat

    @property
    def nsfn(self):
        return self.radial.NSFN()

    @property
    def sfn_lm(self):
        return np.fromiter(
            (item.LM() for item in self.radial.SHAPE_FUNCTIONS.values()),
            dtype=int,
            count=self.nsfn,
        )

    @property
    def sfn(self):
        if not self.nsfn:
            return np.empty((0, self.nr))
        return np.stack(
            [item.VALUES() for item in self.radial.SHAPE_FUNCTIONS.values()]
        )

    @property
    def rmt(self):
        return self.voronoi.RMTRED0() * self._output.alat

    @property
    def rmtfill(self):
        return self.voronoi.RMTFILL()

    @property
    def vol(self):
        return self.voronoi.VOL()

    @property
    def nface(self):
        return self.voronoi.NFACE()

    @property
    def faces(self):
        return tuple(self.voronoi.FACES.values())

    @property
    def face_vertices(self):
        """Vertices of the Voronoi faces, scaled consistently with ``rmesh``."""

        return tuple(
            face.VERTICES() * self._output.alat for face in self.faces
        )

    @property
    def panel_boundaries(self):
        """Radii at the ends of the radial integration panels."""

        indexes = np.cumsum(self.jrcut, dtype=int) - 1
        if len(indexes) and indexes[-1] >= self.nr:
            raise ValueError("Radial panel sizes exceed the radial mesh")
        return self.rmesh[indexes]

    @staticmethod
    def _index2lm(index):
        index = int(index)
        l = math.isqrt(index - 1)
        return l, index - l * l - l - 1

    @staticmethod
    def _real_sph_harm(m, l, phi, theta):
        try:
            from scipy.special import sph_harm_y

            def sph_harmonic(order):
                return sph_harm_y(l, order, theta, phi)

        except ImportError:
            from scipy.special import sph_harm

            def sph_harmonic(order):
                return sph_harm(order, l, phi, theta)

        if m == 0:
            return np.real(sph_harmonic(0))
        if m > 0:
            return np.real(
                (sph_harmonic(-m) + (-1) ** m * sph_harmonic(m))
                / np.sqrt(2)
            )
        return np.real(
            1j
            * (sph_harmonic(-m) - (-1) ** m * sph_harmonic(m))
            / np.sqrt(2)
        )

    def _shape_function(self, x, y, z):
        radius = np.sqrt(x**2 + y**2 + z**2)
        with np.errstate(invalid="ignore", divide="ignore"):
            theta = np.where(
                radius > 0,
                np.arccos(np.clip(z / radius, -1, 1)),
                0.0,
            )
        phi = np.mod(np.arctan2(y, x), 2 * np.pi)

        out = np.zeros_like(radius, dtype=float)
        for values, lm_index in zip(self.sfn, self.sfn_lm):
            l, m = self._index2lm(lm_index)
            radial = np.interp(
                radius,
                self.rmesh,
                values,
                left=values[0],
                right=0.0,
            )
            out += radial * self._real_sph_harm(m, l, phi, theta)
        return out

    def to_3d_grid(self, n=80):
        """Rebuild the shape function on a cubic grid."""

        rmax = self.rmesh[-1] * 1.02
        linear = np.linspace(-rmax, rmax, n)
        x, y, z = np.meshgrid(linear, linear, linear, indexing="ij")
        return linear, self._shape_function(x, y, z)

    def to_2d_grid(self, plane="xy", n=200, offset=0.0):
        """Rebuild a planar section of the shape function."""

        plane = plane.lower()
        if plane not in ("xy", "xz", "yz"):
            raise ValueError("plane has to be one of: xy, xz, yz")
        rmax = self.rmesh[-1] * 1.02
        linear = np.linspace(-rmax, rmax, n)
        first, second = np.meshgrid(linear, linear, indexing="xy")
        fixed = np.full_like(first, offset)
        coordinates = {
            "xy": (first, second, fixed),
            "xz": (first, fixed, second),
            "yz": (fixed, first, second),
        }[plane]
        return linear, self._shape_function(*coordinates)

    def isosurface(self, level=0.5, n=64):
        """Return vertices and triangular faces of a reconstructed isosurface."""

        return self.isosurfaces(levels=(level,), n=n)[0][1:]

    def isosurfaces(self, levels=5, n=64):
        """Return reconstructed isosurfaces as ``(level, vertices, faces)``."""

        try:
            from skimage.measure import marching_cubes
        except ImportError as exc:
            raise ImportError(
                "3D SFN isosurface plotting requires scikit-image; "
                "install ase2sprkkr[plotting] or use plot_slice()"
            ) from exc

        linear, values = self.to_3d_grid(n=n)
        minimum, maximum = np.min(values), np.max(values)
        if np.isscalar(levels):
            count = int(levels)
            if count < 1 or count != levels:
                raise ValueError("levels has to be a positive integer")
            # A truncated spherical-harmonic expansion of the nominally
            # binary shape function may ring below zero or above one.  Those
            # overshoots are not useful automatic contour levels.
            lower = max(float(minimum), 0.0)
            upper = min(float(maximum), 1.0)
            if lower >= upper:
                lower, upper = float(minimum), float(maximum)
            levels = np.linspace(lower, upper, count + 2)[1:-1]
        else:
            levels = np.asarray(tuple(levels), dtype=float)
            if levels.ndim != 1 or not len(levels):
                raise ValueError("levels has to contain at least one value")

        spacing = float(linear[1] - linear[0])
        output = []
        for level in levels:
            if not minimum < level < maximum:
                raise ValueError(
                    f"Isosurface level {level} is outside the open "
                    f"shape-function range ({minimum}, {maximum})"
                )
            vertices, faces, _, _ = marching_cubes(
                values, level=level, spacing=(spacing,) * 3
            )
            vertices += linear[0]
            output.append((float(level), vertices, faces))
        return output

    def plot(self, what="mesh", **kwargs):
        """Plot the Voronoi mesh, a shape-function section or radial data."""

        methods = {
            "mesh": self.plot_mesh,
            "shape": self.plot_shape,
            "slice": self.plot_slice,
            "radial": self.plot_radial,
        }
        try:
            method = methods[what.lower()]
        except (AttributeError, KeyError) as exc:
            raise ValueError(
                "what has to be one of: mesh, shape, slice, radial"
            ) from exc
        return method(**kwargs)

    def plot_mesh(
        self,
        axis=None,
        *,
        show_rmt=True,
        show_panels=False,
        facecolor="tab:blue",
        alpha=0.22,
        edgecolor="black",
        filename=None,
        show=None,
        dpi=600,
        latex=None,
        figsize=(6, 5),
    ):
        """Plot the Voronoi polyhedron and optional radial boundaries."""

        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        with single_plot(
            axis=axis,
            projection="3d",
            filename=filename,
            show=show,
            dpi=dpi,
            latex=latex,
            figsize=figsize,
        ) as (axis, _):
            vertices = self.face_vertices
            collection = Poly3DCollection(
                vertices,
                facecolor=facecolor,
                edgecolor=edgecolor,
                alpha=alpha,
            )
            axis.add_collection3d(collection)

            sphere_radii = []
            if show_rmt:
                sphere_radii.append((self.rmt, 0.18, "tab:orange"))
            if show_panels:
                sphere_radii.extend(
                    (radius, 0.08, "0.45")
                    for radius in self.panel_boundaries[:-1]
                    if not np.isclose(radius, self.rmt)
                )
            u = np.linspace(0, 2 * np.pi, 32)
            v = np.linspace(0, np.pi, 16)
            for radius, sphere_alpha, color in sphere_radii:
                x = radius * np.outer(np.cos(u), np.sin(v))
                y = radius * np.outer(np.sin(u), np.sin(v))
                z = radius * np.outer(np.ones_like(u), np.cos(v))
                axis.plot_surface(
                    x,
                    y,
                    z,
                    color=color,
                    alpha=sphere_alpha,
                    linewidth=0,
                )

            all_vertices = np.concatenate(vertices)
            extent = np.max(np.abs(all_vertices)) * 1.08
            for setter in (axis.set_xlim, axis.set_ylim, axis.set_zlim):
                setter(-extent, extent)
            axis.set_box_aspect((1, 1, 1))
            axis.set_xlabel("x")
            axis.set_ylabel("y")
            axis.set_zlabel("z")
            axis.set_title(f"Shape-function mesh IM={self.idx}")
        return axis

    def plot_shape(
        self,
        axis=None,
        *,
        n=None,
        levels=1,
        level=None,
        facecolor=None,
        cmap="viridis",
        alpha=0.25,
        edgecolor="none",
        filename=None,
        show=None,
        dpi=600,
        latex=None,
        figsize=(6, 5),
    ):
        """Plot isosurfaces of the reconstructed shape function."""

        with single_plot(
            axis=axis,
            projection="3d",
            filename=filename,
            show=show,
            dpi=dpi,
            latex=latex,
            figsize=figsize,
        ) as (axis, _):
            n = 64 if n is None else n
            surfaces = self.isosurfaces(
                levels=(level,) if level is not None else levels, n=n
            )

            from mpl_toolkits.mplot3d.art3d import Poly3DCollection
            import matplotlib.pyplot as plt

            colors = (
                [facecolor] * len(surfaces)
                if facecolor is not None
                else plt.get_cmap(cmap)(np.linspace(0.1, 0.9, len(surfaces)))
            )
            for color, (_, vertices, faces) in zip(colors, surfaces):
                collection_options = {
                    "facecolors": color,
                    "alpha": alpha,
                }
                if edgecolor not in (None, "none"):
                    collection_options["edgecolors"] = edgecolor
                axis.add_collection3d(
                    Poly3DCollection(vertices[faces], **collection_options)
                )
            all_vertices = np.concatenate(
                [vertices for _, vertices, _ in surfaces]
            )
            extent = np.max(np.abs(all_vertices)) * 1.05
            for setter in (axis.set_xlim, axis.set_ylim, axis.set_zlim):
                setter(-extent, extent)
            axis.set_box_aspect((1, 1, 1))
            axis.set_xlabel("x")
            axis.set_ylabel("y")
            axis.set_zlabel("z")
            axis.set_title(
                f"Shape-function isosurfaces IM={self.idx} "
                f"({len(surfaces)} levels)"
            )
        return axis

    def plot_slice(
        self,
        axis=None,
        *,
        plane="xy",
        n=200,
        offset=0.0,
        level=0.5,
        levels=31,
        cmap="viridis",
        colorbar=True,
        show_boundary=True,
        filename=None,
        show=None,
        dpi=600,
        latex=None,
        figsize=(6, 5),
    ):
        """Plot a planar section through the reconstructed shape function."""

        with single_plot(
            axis=axis,
            filename=filename,
            show=show,
            dpi=dpi,
            latex=latex,
            figsize=figsize,
        ) as (axis, _):
            linear, values = self.to_2d_grid(plane=plane, n=n, offset=offset)
            contour = axis.contourf(linear, linear, values, levels=levels, cmap=cmap)
            if show_boundary and np.min(values) <= level <= np.max(values):
                axis.contour(
                    linear,
                    linear,
                    values,
                    levels=[level],
                    colors="black",
                    linewidths=1,
                )
            if colorbar:
                axis.figure.colorbar(contour, ax=axis, label="shape function")
            labels = tuple(plane.lower())
            axis.set_xlabel(labels[0])
            axis.set_ylabel(labels[1])
            axis.set_aspect("equal")
            axis.set_title(f"Shape function IM={self.idx}, {plane.lower()} section")
        return axis

    def plot_radial(
        self,
        axis=None,
        *,
        show_panels=True,
        legend=True,
        filename=None,
        show=None,
        dpi=600,
        latex=None,
        figsize=(7, 4),
    ):
        """Plot the radial spherical-harmonic coefficients."""

        with single_plot(
            axis=axis,
            filename=filename,
            show=show,
            dpi=dpi,
            latex=latex,
            figsize=figsize,
        ) as (axis, _):
            for values, lm_index in zip(self.sfn, self.sfn_lm):
                l, m = self._index2lm(lm_index)
                axis.plot(self.rmesh, values, label=rf"$l={l}, m={m}$")
            if show_panels:
                for radius in self.panel_boundaries[:-1]:
                    axis.axvline(radius, color="0.6", linewidth=0.7, linestyle="--")
            axis.set_xlabel("r")
            axis.set_ylabel(r"$S_{lm}(r)$")
            axis.set_title(f"Radial shape functions IM={self.idx}")
            axis.grid(alpha=0.25)
            if legend:
                axis.legend(fontsize="small", ncols=min(3, max(1, self.nsfn)))
        return axis

    def __repr__(self):
        return (
            f"ShapeFunctionMesh(idx={self.idx}, npan={self.npan}, "
            f"nr={self.nr}, nsfn={self.nsfn}, nface={self.nface})"
        )


class _MeshesView(Sequence):
    def __init__(self, output):
        self.output = output

    def __len__(self):
        return self.output.NM()

    def __getitem__(self, index):
        if isinstance(index, slice):
            return tuple(ShapeFunctionMesh(self.output, i) for i in range(*index.indices(len(self))))
        if index < 0:
            index += len(self)
        if index < 0 or index >= len(self):
            raise IndexError(index)
        return ShapeFunctionMesh(self.output, index)


class SFNOutputFile(OutputFile):
    """Grammar-backed representation of an SPR-KKR shape-function file."""

    additional_actions = (
        "plot_cell",
        "plot_periodic",
        "plot_shape",
        "plot_slice",
        "plot_radial",
    )
    plot_parameters = {
        "what",
        "mesh",
        "plane",
        "n",
        "offset",
        "level",
        "levels",
        "cmap",
        "conventional",
        "separate_plots",
    }

    @property
    def alat(self):
        """Length scale in Angstrom, taken from the associated potential."""

        from ase.units import Bohr

        return self.potential.LATTICE.ALAT() * Bohr

    @property
    def nm(self):
        """Return the number of shape-function meshes."""

        return self.NM()

    @property
    def meshes(self):
        """Provide indexed access to shape-function mesh views."""

        return _MeshesView(self)

    def mesh_for_idx(self, index):
        """Return the mesh with the given SPR-KKR index, if present."""

        return next((mesh for mesh in self.meshes if mesh.idx == index), None)

    def _selected_meshes(self, mesh):
        """Resolve a mesh selector to a non-empty tuple of mesh views."""

        if mesh is None:
            selected = tuple(self.meshes)
        elif isinstance(mesh, ShapeFunctionMesh):
            selected = (mesh,)
        elif isinstance(mesh, slice):
            selected = self.meshes[mesh]
        elif isinstance(mesh, (int, np.integer)):
            selected = (self.meshes[int(mesh)],)
        else:
            try:
                selected = tuple(self.meshes[int(index)] for index in mesh)
            except (TypeError, ValueError) as exc:
                raise TypeError(
                    "mesh has to be an index, slice, ShapeFunctionMesh "
                    "or iterable of indexes"
                ) from exc
        if not selected:
            raise ValueError("No shape-function mesh selected")
        return selected

    def _cell_context(self, selected, conventional):
        """Return atomic, mesh-index and cell data used by cell plots."""

        try:
            potential = self.potential
            atoms = potential.atoms
        except (OSError, ValueError) as exc:
            raise ValueError(
                "A potential is required to plot SFN data in the unit cell. "
                "Pass potential=... to OutputFile.from_file(), or place the "
                "matching .pot or .pot_new file beside the SFN file."
            ) from exc
        if atoms is None:
            raise ValueError("The SFN potential does not contain an atomic structure")

        rows = potential.OCCUPATION.DATA()
        if len(rows) != len(atoms):
            raise ValueError(
                "The number of OCCUPATION rows does not match the number of atoms"
            )
        all_meshes = {item.idx for item in self.meshes}
        invalid = set(np.asarray(rows["IMQ"], dtype=int)) - all_meshes
        if invalid:
            raise ValueError(
                "The potential refers to undefined SFN mesh(es): "
                + ", ".join(map(str, sorted(invalid)))
            )

        mesh_indices = np.asarray(rows["IMQ"], dtype=int)
        if conventional:
            from ase import Atoms
            from ase.build import make_supercell

            displayed_atoms = Atoms(
                numbers=atoms.get_atomic_numbers(),
                positions=atoms.positions,
                cell=atoms.cell,
                pbc=True,
            )
            displayed_atoms.set_array("sfn_mesh", mesh_indices)
            transformation = np.asarray(
                atoms.cell.get_bravais_lattice().conventional_cellmap,
                dtype=int,
            )
            atoms = make_supercell(displayed_atoms, transformation, wrap=True)
            mesh_indices = atoms.get_array("sfn_mesh")

        cell = np.asarray(atoms.cell, dtype=float)
        if abs(np.linalg.det(cell)) < 1e-12:
            raise ValueError("The potential does not define a three-dimensional cell")
        return (
            atoms,
            mesh_indices,
            {item.idx: item for item in selected},
            cell,
            np.linalg.inv(cell),
        )

    @staticmethod
    def _clip_polygon_to_cell(vertices, cell, inverse_cell):
        """Clip a Cartesian polygon to the unit-cell parallelepiped."""

        polygon = np.asarray(vertices, dtype=float) @ inverse_cell
        if np.any(polygon.max(axis=0) < -1e-10) or np.any(
            polygon.min(axis=0) > 1 + 1e-10
        ):
            return None
        if np.all(polygon >= -1e-10) and np.all(polygon <= 1 + 1e-10):
            return np.clip(polygon, 0.0, 1.0) @ cell
        for dimension in range(3):
            for boundary, keep_greater in ((0.0, True), (1.0, False)):
                if not len(polygon):
                    return None
                clipped = []
                previous = polygon[-1]
                previous_inside = (
                    previous[dimension] >= boundary - 1e-10
                    if keep_greater
                    else previous[dimension] <= boundary + 1e-10
                )
                for current in polygon:
                    current_inside = (
                        current[dimension] >= boundary - 1e-10
                        if keep_greater
                        else current[dimension] <= boundary + 1e-10
                    )
                    if current_inside != previous_inside:
                        delta = current[dimension] - previous[dimension]
                        if abs(delta) > 1e-15:
                            fraction = (boundary - previous[dimension]) / delta
                            clipped.append(previous + fraction * (current - previous))
                    if current_inside:
                        clipped.append(current)
                    previous = current
                    previous_inside = current_inside
                polygon = np.asarray(clipped)
        if len(polygon) < 3:
            return None
        return np.clip(polygon, 0.0, 1.0) @ cell

    @staticmethod
    def _cell_corners(cell):
        """Return the eight Cartesian corners of a unit cell."""

        return np.asarray(
            [
                i * cell[0] + j * cell[1] + k * cell[2]
                for i in (0, 1)
                for j in (0, 1)
                for k in (0, 1)
            ]
        )

    @staticmethod
    def _draw_cell(axis, cell, color):
        """Draw unit-cell edges and return their Cartesian corners."""

        corners = SFNOutputFile._cell_corners(cell)
        for start in range(8):
            for bit in (1, 2, 4):
                end = start ^ bit
                if start < end:
                    axis.plot(
                        *corners[[start, end]].T,
                        color=color,
                        linewidth=1.2,
                    )
        return corners

    @staticmethod
    def _draw_polygon_edges(axis, polygons, color, linewidth):
        """Draw polygons of unequal sizes without confusing Matplotlib autoscaling."""

        from mpl_toolkits.mplot3d.art3d import Line3DCollection

        grouped = {}
        for polygon in polygons:
            line = np.vstack((polygon, polygon[0]))
            grouped.setdefault(len(line), []).append(line)
        for lines in grouped.values():
            axis.add_collection3d(
                Line3DCollection(lines, colors=color, linewidths=linewidth)
            )

    @staticmethod
    def _periodic_centers(position, vertices, cell, inverse_cell):
        """Yield translations whose geometry intersects the displayed cell."""

        center = np.asarray(position) @ inverse_cell
        relative = np.asarray(vertices) @ inverse_cell
        lower = np.ceil(-relative.max(axis=0) - center - 1e-10).astype(int)
        upper = np.floor(1 - relative.min(axis=0) - center + 1e-10).astype(int)
        for i in range(lower[0], upper[0] + 1):
            for j in range(lower[1], upper[1] + 1):
                for k in range(lower[2], upper[2] + 1):
                    yield np.asarray(position) + np.asarray((i, j, k)) @ cell

    @staticmethod
    def _periodic_mesh_colors(meshes, facecolor, colormap="hsv"):
        """Assign a distinct plotting color to every selected mesh."""

        if facecolor is not None:
            return {index: facecolor for index in meshes}

        import matplotlib.pyplot as plt

        colormap = plt.get_cmap(colormap)
        indexes = sorted(meshes)
        return {
            index: colormap((position / len(indexes) + 0.04) % 1.0)
            for position, index in enumerate(indexes)
        }

    def _plot_cell_on_axis(
        self,
        axis,
        selected,
        *,
        show_atoms=True,
        conventional=True,
        facecolor=None,
        alpha=0.35,
        edgecolor=None,
        cellcolor="black",
        legend=True,
        colormap="hsv",
    ):
        """Plot periodically tiled and cell-clipped SFN geometry."""

        atoms, mesh_indices, meshes, cell, inverse_cell = self._cell_context(
            selected, conventional
        )
        positions = atoms.get_scaled_positions(wrap=True) @ cell
        colors = self._periodic_mesh_colors(meshes, facecolor, colormap)

        geometries = {}
        for mesh in selected:
            polygons = mesh.face_vertices
            geometries[mesh.idx] = (np.concatenate(polygons), polygons)

        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        polygons = {index: [] for index in meshes}
        for position, mesh_index in zip(positions, mesh_indices):
            mesh_index = int(mesh_index)
            if mesh_index not in meshes:
                continue
            extent_vertices, source_polygons = geometries[mesh_index]
            for center in self._periodic_centers(
                position, extent_vertices, cell, inverse_cell
            ):
                for source in source_polygons:
                    polygon = self._clip_polygon_to_cell(
                        source + center, cell, inverse_cell
                    )
                    if polygon is not None:
                        polygons[mesh_index].append(polygon)
        if not any(polygons.values()):
            raise ValueError("The selected SFN geometry does not intersect the cell")
        for mesh_index, mesh_polygons in polygons.items():
            if not mesh_polygons:
                continue
            color = colors[mesh_index]
            axis.add_collection3d(
                Poly3DCollection(
                    mesh_polygons,
                    facecolor=color,
                    edgecolor=edgecolor or color,
                    alpha=alpha,
                    linewidths=0.6,
                )
            )

        corners = self._draw_cell(axis, cell, cellcolor)
        if show_atoms:
            shown_positions = []
            shown_colors = []
            for position, mesh_index in zip(positions, mesh_indices):
                mesh_index = int(mesh_index)
                if mesh_index not in colors:
                    continue
                fractional = position @ inverse_cell
                lower = np.ceil(-fractional - 1e-10).astype(int)
                upper = np.floor(1 - fractional + 1e-10).astype(int)
                for i in range(lower[0], upper[0] + 1):
                    for j in range(lower[1], upper[1] + 1):
                        for k in range(lower[2], upper[2] + 1):
                            shown_positions.append(
                                position + np.asarray((i, j, k)) @ cell
                            )
                            shown_colors.append(colors[mesh_index])
            shown_positions = np.asarray(shown_positions)
            axis.scatter(
                *shown_positions.T,
                color=shown_colors,
                edgecolor="black",
                linewidth=0.5,
                s=70,
                depthshade=True,
            )

        if legend and facecolor is None:
            from ase.data import chemical_symbols
            from matplotlib.patches import Patch

            symbols = {}
            for number, mesh_index in zip(atoms.numbers, mesh_indices):
                symbols.setdefault(int(mesh_index), set()).add(
                    chemical_symbols[number]
                )
            handles = [
                Patch(
                    facecolor=colors[index],
                    edgecolor=colors[index],
                    alpha=max(alpha, 0.5),
                    label=(
                        f"IMQ {index}: "
                        + ", ".join(sorted(symbols.get(index, ())))
                    ).rstrip(": "),
                )
                for index in meshes
            ]
            axis.legend(handles=handles, loc="upper left", fontsize="small")

        minimum = corners.min(axis=0)
        maximum = corners.max(axis=0)
        center = (minimum + maximum) / 2
        extent = max(maximum - minimum) / 2
        if extent == 0:
            extent = 0.5
        for setter, value in zip(
            (axis.set_xlim, axis.set_ylim, axis.set_zlim), center
        ):
            setter(value - extent, value + extent)
        axis.set_box_aspect((1, 1, 1))
        axis.set_proj_type("ortho")
        axis.view_init(elev=22, azim=-55)
        axis.grid(False)
        for coordinate_axis in (axis.xaxis, axis.yaxis, axis.zaxis):
            coordinate_axis.pane.fill = False
            coordinate_axis.pane.set_edgecolor("none")
        axis.set_xlabel("x [Angstrom]")
        axis.set_ylabel("y [Angstrom]")
        axis.set_zlabel("z [Angstrom]")
        kind = "conventional" if conventional else "primitive"
        axis.set_title(f"SFN meshes in the {kind} cell")

    def _plot_representative_on_axis(
        self,
        axis,
        mesh,
        cell,
        *,
        shape,
        show_atom=True,
        show_mesh=True,
        n=64,
        levels=1,
        level=None,
        facecolor="tab:blue",
        cmap="viridis",
        alpha=0.22,
        edgecolor="tab:blue",
        cellcolor="black",
        conventional=True,
    ):
        """Draw one unique SFN mesh without overlapping periodic copies."""

        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        if shape:
            import matplotlib.pyplot as plt

            surfaces = mesh.isosurfaces(
                levels=(level,) if level is not None else levels, n=n
            )
            colors = (
                [facecolor] * len(surfaces)
                if facecolor is not None
                else plt.get_cmap(cmap)(np.linspace(0.1, 0.9, len(surfaces)))
            )
            for color, (_, vertices, faces) in zip(colors, surfaces):
                axis.add_collection3d(
                    Poly3DCollection(
                        vertices[faces],
                        facecolors=color,
                        alpha=alpha,
                    )
                )
            vertices = np.concatenate(
                [vertices for _, vertices, _ in surfaces]
            )
        else:
            vertices = np.concatenate(mesh.face_vertices)
            polygons = mesh.face_vertices
            collection_options = {
                "facecolors": facecolor,
                "linewidths": 0.7,
                "alpha": alpha,
            }
            if edgecolor not in (None, "none"):
                collection_options["edgecolors"] = edgecolor
            axis.add_collection3d(
                Poly3DCollection(polygons, **collection_options)
            )
        if shape and show_mesh:
            self._draw_polygon_edges(axis, mesh.face_vertices, edgecolor, 0.8)

        origin = -0.5 * cell.sum(axis=0)
        corners = self._draw_cell(axis, cell, cellcolor) + origin
        # _draw_cell draws an origin-based cell, so translate its lines as well.
        for line in axis.lines[-12:]:
            x, y, z = line.get_data_3d()
            translated = np.column_stack((x, y, z)) + origin
            line.set_data_3d(*translated.T)

        if show_atom:
            axis.scatter(
                0,
                0,
                0,
                color="0.65",
                edgecolor="black",
                linewidth=0.6,
                s=85,
                depthshade=True,
            )

        points = np.concatenate((corners, vertices))
        minimum = points.min(axis=0)
        maximum = points.max(axis=0)
        center = (minimum + maximum) / 2
        extent = max(maximum - minimum) / 2 * 1.05
        for setter, value in zip(
            (axis.set_xlim, axis.set_ylim, axis.set_zlim), center
        ):
            setter(value - extent, value + extent)
        axis.set_box_aspect((1, 1, 1))
        axis.set_proj_type("ortho")
        axis.view_init(elev=22, azim=-55)
        axis.grid(False)
        for coordinate_axis in (axis.xaxis, axis.yaxis, axis.zaxis):
            coordinate_axis.pane.fill = False
            coordinate_axis.pane.set_edgecolor("none")
        axis.set_xlabel("x [Angstrom]")
        axis.set_ylabel("y [Angstrom]")
        axis.set_zlabel("z [Angstrom]")
        kind = "conventional" if conventional else "primitive"
        subject = "isosurface" if shape else "mesh"
        if shape:
            level_text = f", {len(surfaces)} levels"
        else:
            level_text = ""
        axis.set_title(
            f"SFN {subject} IM={mesh.idx}{level_text}; centered {kind} cell"
        )

    def plot(self, what="mesh", **kwargs):
        """Dispatch to one of the explicit SFN visualization methods."""

        methods = {
            "mesh": self.plot_mesh,
            "cell": self.plot_cell,
            "periodic": self.plot_periodic,
            "shape": self.plot_shape,
            "slice": self.plot_slice,
            "radial": self.plot_radial,
        }
        try:
            method = methods[what.lower()]
        except (AttributeError, KeyError) as exc:
            raise ValueError(
                "what has to be one of: mesh, cell, periodic, shape, slice, radial"
            ) from exc
        return method(**kwargs)

    def plot_mesh(
        self,
        *,
        mesh=None,
        layout=None,
        separate_plots=False,
        figsize=None,
        filename=None,
        show=None,
        dpi=600,
        latex=None,
        show_rmt=True,
        show_panels=False,
        facecolor="tab:blue",
        alpha=0.22,
        edgecolor="black",
    ):
        """Plot one or more Voronoi meshes."""

        selected = self._selected_meshes(mesh)
        if figsize is None:
            figsize = lambda size: (6 * size[1], 5 * size[0])
        with Multiplot(
            layout=layout,
            number_of_plots=len(selected),
            projection="3d",
            separate_plots=separate_plots,
            figsize=figsize,
            filename=filename,
            show=show,
            dpi=dpi,
            latex=latex,
        ) as multiplot:
            for item in selected:
                multiplot.plot(
                    item,
                    name=f"IM{item.idx}",
                    plot_function=item.plot_mesh,
                    show_rmt=show_rmt,
                    show_panels=show_panels,
                    facecolor=facecolor,
                    alpha=alpha,
                    edgecolor=edgecolor,
                )
        return multiplot.figure

    def plot_cell(
        self,
        *,
        mesh=None,
        layout=None,
        separate_plots=False,
        figsize=None,
        filename=None,
        show=None,
        dpi=600,
        latex=None,
        show_atoms=True,
        facecolor="tab:blue",
        alpha=0.16,
        edgecolor="tab:blue",
        conventional=True,
        cellcolor="black",
    ):
        """Plot each unique mesh in a centered crystallographic cell."""

        selected = self._selected_meshes(mesh)
        _, _, _, cell, _ = self._cell_context(selected, conventional)
        if figsize is None:
            figsize = lambda size: (6 * size[1], 5 * size[0])
        with Multiplot(
            layout=layout,
            number_of_plots=len(selected),
            projection="3d",
            separate_plots=separate_plots,
            figsize=figsize,
            filename=filename,
            show=show,
            dpi=dpi,
            latex=latex,
        ) as multiplot:
            for item in selected:

                def plot_function(axis, item=item):
                    self._plot_representative_on_axis(
                        axis,
                        item,
                        cell,
                        shape=False,
                        show_atom=show_atoms,
                        facecolor=facecolor,
                        alpha=alpha,
                        edgecolor=edgecolor,
                        cellcolor=cellcolor,
                        conventional=conventional,
                    )

                multiplot.plot(
                    item,
                    name=f"IM{item.idx}",
                    plot_function=plot_function,
                )
        return multiplot.figure

    def plot_periodic(
        self,
        *,
        mesh=None,
        layout=None,
        separate_plots=False,
        figsize=None,
        filename=None,
        show=None,
        dpi=600,
        latex=None,
        show_atoms=True,
        facecolor=None,
        alpha=0.16,
        edgecolor=None,
        conventional=True,
        cellcolor="black",
        legend=True,
        colormap="hsv",
    ):
        """Plot the periodically tiled meshes clipped to one cell."""

        selected = self._selected_meshes(mesh)
        with Multiplot(
            layout=layout,
            number_of_plots=1,
            projection="3d",
            separate_plots=separate_plots,
            figsize=figsize or (7, 6),
            filename=filename,
            show=show,
            dpi=dpi,
            latex=latex,
        ) as multiplot:

            def plot_function(axis):
                self._plot_cell_on_axis(
                    axis,
                    selected,
                    show_atoms=show_atoms,
                    conventional=conventional,
                    facecolor=facecolor,
                    alpha=alpha,
                    edgecolor=edgecolor,
                    cellcolor=cellcolor,
                    legend=legend,
                    colormap=colormap,
                )

            multiplot.plot(self, name="cell", plot_function=plot_function)
        return multiplot.figure

    def plot_shape(
        self,
        *,
        mesh=None,
        layout=None,
        separate_plots=False,
        figsize=None,
        filename=None,
        show=None,
        dpi=600,
        latex=None,
        n=None,
        levels=1,
        level=None,
        facecolor=None,
        cmap="viridis",
        alpha=0.25,
        edgecolor="none",
        show_atoms=True,
        show_mesh=True,
        conventional=True,
        cellcolor="black",
    ):
        """Plot 3D shape-function isosurfaces for each unique mesh."""

        selected = self._selected_meshes(mesh)
        _, _, _, cell, _ = self._cell_context(selected, conventional)
        if figsize is None:
            figsize = lambda size: (6 * size[1], 5 * size[0])
        with Multiplot(
            layout=layout,
            number_of_plots=len(selected),
            projection="3d",
            separate_plots=separate_plots,
            figsize=figsize,
            filename=filename,
            show=show,
            dpi=dpi,
            latex=latex,
        ) as multiplot:
            for item in selected:

                def plot_function(axis, item=item):
                    self._plot_representative_on_axis(
                        axis,
                        item,
                        cell,
                        shape=True,
                        show_atom=show_atoms,
                        show_mesh=show_mesh,
                        n=64 if n is None else n,
                        levels=levels,
                        level=level,
                        facecolor=facecolor,
                        cmap=cmap,
                        alpha=alpha,
                        edgecolor=(
                            "black" if edgecolor == "none" else edgecolor
                        ),
                        cellcolor=cellcolor,
                        conventional=conventional,
                    )

                multiplot.plot(
                    item,
                    name=f"IM{item.idx}",
                    plot_function=plot_function,
                )
        return multiplot.figure

    def plot_slice(
        self,
        *,
        mesh=None,
        layout=None,
        separate_plots=False,
        figsize=None,
        filename=None,
        show=None,
        dpi=600,
        latex=None,
        plane="xy",
        n=200,
        offset=0.0,
        level=0.5,
        levels=31,
        cmap="viridis",
        colorbar=True,
        show_boundary=True,
    ):
        """Plot planar sections through one or more shape functions."""

        selected = self._selected_meshes(mesh)
        if figsize is None:
            figsize = lambda size: (6 * size[1], 5 * size[0])
        with Multiplot(
            layout=layout,
            number_of_plots=len(selected),
            separate_plots=separate_plots,
            figsize=figsize,
            filename=filename,
            show=show,
            dpi=dpi,
            latex=latex,
        ) as multiplot:
            for item in selected:
                multiplot.plot(
                    item,
                    name=f"IM{item.idx}",
                    plot_function=item.plot_slice,
                    plane=plane,
                    n=n,
                    offset=offset,
                    level=level,
                    levels=levels,
                    cmap=cmap,
                    colorbar=colorbar,
                    show_boundary=show_boundary,
                )
        return multiplot.figure

    def plot_radial(
        self,
        *,
        mesh=None,
        layout=None,
        separate_plots=False,
        figsize=None,
        filename=None,
        show=None,
        dpi=600,
        latex=None,
        show_panels=True,
        legend=True,
    ):
        """Plot radial coefficients of one or more shape functions."""

        selected = self._selected_meshes(mesh)
        if figsize is None:
            figsize = lambda size: (6 * size[1], 4 * size[0])
        with Multiplot(
            layout=layout,
            number_of_plots=len(selected),
            separate_plots=separate_plots,
            figsize=figsize,
            filename=filename,
            show=show,
            dpi=dpi,
            latex=latex,
        ) as multiplot:
            for item in selected:
                multiplot.plot(
                    item,
                    name=f"IM{item.idx}",
                    plot_function=item.plot_radial,
                    show_panels=show_panels,
                    legend=legend,
                )
        return multiplot.figure

    def __repr__(self):
        try:
            alat = f"{self.alat:.6g}"
        except (AttributeError, OSError, ValueError):
            alat = "unknown"
        return f"SFNOutputFile(nm={self.nm}, alat={alat})"


class SFNDefinition(OutputFileDefinition):
    result_class = SFNOutputFile


def create_definition():
    V = _SFNValueDefinition
    S = _SFNSectionDefinition
    GV = GeneratedValueDefinition

    header = V(
        "HEADER",
        NumpyArray(
            delimiter=5,
            shape=(3,),
            items_per_line=3,
            item_format="%5d",
            written_delimiter="",
            dtype=int,
        ),
        name_in_grammar=False,
    )

    radial_mesh = S(
        "RADIAL_MESHES",
        [
            V(
                "MESH_COUNTS",
                NumpyArray(
                    delimiter=5,
                    shape=(2,),
                    items_per_line=2,
                    item_format="%5d",
                    written_delimiter="",
                    dtype=int,
                ),
                name_in_grammar=False,
            ),
            GV("NPAN", lambda section: int(section.MESH_COUNTS()[0])),
            GV("NR", lambda section: int(section.MESH_COUNTS()[1])),
            V(
                "JRCUT",
                NumpyArray(
                    delimiter=5,
                    shape=_ArrayShape.as_dependent_value("MESH_COUNTS", 0),
                    items_per_line=16,
                    item_format="%5d",
                    written_delimiter="",
                    dtype=int,
                ),
                name_in_grammar=False,
            ),
            V(
                "RMESH",
                NumpyArray(
                    delimiter=22,
                    shape=_ArrayShape.as_dependent_value("MESH_COUNTS", 1, 2),
                    items_per_line=4,
                    item_format="%22.14E",
                    written_delimiter="",
                    dtype=float,
                ),
                name_in_grammar=False,
            ),
            V("NSFN", Integer(format="5d"), name_in_grammar=False),
            S(
                "SHAPE_FUNCTIONS",
                [
                    V("LM", Integer(format="5d"), name_in_grammar=False),
                    V(
                        "VALUES",
                        NumpyArray(
                            delimiter=22,
                            shape=_ArrayShape.as_dependent_value("..MESH_COUNTS", 1),
                            items_per_line=4,
                            item_format="%22.14E",
                            written_delimiter="",
                            dtype=float,
                        ),
                        name_in_grammar=False,
                    ),
                ],
                name_in_grammar=False,
                is_repeated=True,
                repeated_count="NSFN",
                repeated_with_name=False,
                repeated_delimiter=line_end,
            ),
        ],
        name_in_grammar=False,
        is_repeated=True,
        repeated_count=_ArrayItem.as_dependent_value("HEADER", 0),
        repeated_with_name=False,
        repeated_delimiter=line_end,
    )

    face = S(
        "FACES",
        [
            V(
                "FACE_COUNTS",
                NumpyArray(
                    delimiter=10,
                    shape=(2,),
                    items_per_line=2,
                    item_format="%10d",
                    written_delimiter="",
                    dtype=int,
                ),
                written_name="IFC NVERT_FCM(IFC,IM)",
                delimiter=line_end,
            ),
            GV("IFC", lambda section: int(section.FACE_COUNTS()[0])),
            GV("NVERT", lambda section: int(section.FACE_COUNTS()[1])),
            V(
                "VERTICES",
                NumpyArray(
                    delimiter=25,
                    shape=_ArrayShape.as_dependent_value("FACE_COUNTS", 1, 3),
                    items_per_line=3,
                    item_format="%25.17E",
                    written_delimiter="",
                    dtype=float,
                ),
                name_in_grammar=False,
            ),
            V(
                "ACOEF",
                NumpyArray(
                    delimiter=25,
                    shape=(4,),
                    items_per_line=4,
                    item_format="%25.17E",
                    written_delimiter="",
                    dtype=float,
                ),
                written_name="ACOEF_FM",
                delimiter=line_end,
            ),
        ],
        name_in_grammar=False,
        is_repeated=True,
        repeated_count="NFACE",
        repeated_with_name=False,
        repeated_delimiter=line_end,
    )

    voronoi_mesh = S(
        "VORONOI_MESHES",
        [
            *_fortran_record(
                V(
                    "IM",
                    Integer(format="10d"),
                    written_name="IM RMTRED0_M(IM) RMTFILL_M(IM) VOL_M(IM)",
                    delimiter=line_end,
                ),
                V("RMTRED0", Real(format="25.17E"), name_in_grammar=False),
                V("RMTFILL", Real(format="25.17E"), name_in_grammar=False),
                V("VOL", Real(format="25.17E"), name_in_grammar=False),
            ),
            V(
                "NFACE",
                Integer(format="10d"),
                written_name="NFACE_M(IM)",
                delimiter=line_end,
            ),
            face,
        ],
        name_in_grammar=False,
        is_repeated=True,
        repeated_count=_ArrayItem.as_dependent_value("HEADER", 0),
        repeated_with_name=False,
        repeated_delimiter=line_end,
    )

    return SFNDefinition(
        "SFN",
        [
            header,
            GV("NM", lambda section: int(section.HEADER()[0])),
            GV("IFMTSFN", lambda section: int(section.HEADER()[1])),
            GV("NL", lambda section: int(section.HEADER()[2])),
            radial_mesh,
            voronoi_mesh,
        ],
    )


definition = create_definition()


__all__ = ["SFNOutputFile", "ShapeFunctionMesh", "SFNDefinition", "definition"]
