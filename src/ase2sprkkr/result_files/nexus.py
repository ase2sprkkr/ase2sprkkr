"""Export native SPR-KKR ARPES SPC data as NeXus/HDF5.

The SPC file owns sampled intensities and coordinates. InputParameters owns
the *specified* experiment settings. This module does not recompute KKRSPEC
spin, photon, or reciprocal-space transformations.
"""

from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version as distribution_version
from pathlib import Path
import os
import tempfile
import warnings

import numpy as np


# A single declarative map for optional, *original* SPR-KKR settings.
# Values in this map are not inferred effective values.
_INPUT_OPTIONS = {
    "TASK": ("TASK", "INPVER", "STRVER", "IQ_AT_SURF",
             "MILLER_HKL", "CRYS_VECS"),
    "ENERGY": ("EMIN", "EMAX", "NE", "EWORK_EV", "IMV_INI_EV", "IMV_FIN_EV"),
    "SPEC_PH": ("EPHOT", "THETA", "PHI", "POL_P", "NPOL", "ICIRC",
                "IDREH", "ALQ", "DELQ", "IFSP", "THETA_FIX"),
    "SPEC_EL": ("THETA", "PHI", "NT", "NP", "KA",
                "K1", "K2", "K3", "K4", "NK1", "NK2", "NK3", "NK4",
                "SPOL", "PSPIN", "POL_E", "TYP", "ISTR", "POL0", "POL0L",
                "Q1", "Q2", "Q3", "Q4", "BETA1", "BETA2", "ROTAXIS"),
    "SPEC_STR": ("SURF_BAR", "TRANSP_BAR", "N_LAYDBL",
                 "NLAT_G_VEC", "N_LAYER"),
}
_EXPORT_METADATA = {
    "photon_energy": ("SPEC_PH", "EPHOT"),
    "photon_theta": ("SPEC_PH", "THETA"),
    "photon_phi": ("SPEC_PH", "PHI"),
    "photon_polarization": ("SPEC_PH", "POL_P"),
    "work_function": ("ENERGY", "EWORK_EV"),
    "surface_hkl": ("TASK", "MILLER_HKL"),
}

# Units of original input parameters are added only where their semantics
# are established by the parameter names/definitions. Do not infer effective
# photon or electron geometry from the input.
_INPUT_UNITS = {
    ("ENERGY", "EWORK_EV"): "eV",
    ("ENERGY", "IMV_INI_EV"): "eV",
    ("ENERGY", "IMV_FIN_EV"): "eV",
    ("SPEC_PH", "EPHOT"): "eV",
    ("SPEC_PH", "THETA"): "degree",
    ("SPEC_PH", "PHI"): "degree",
    ("SPEC_EL", "THETA"): "degree",
    ("SPEC_EL", "PHI"): "degree",
}
# Native SPC header fields that the parser exposes directly. The original
# values are preserved, without guessing the units of EFERMI or applying
# physical conversions.
_SPC_HEADER = {
    "keyword": "KEYWORD",
    "title": "TITLE",
    "system": "SYSTEM",
    "number_of_sites": "NQ_EFF",
    "number_of_atom_types": "NT_EFF",
    "number_of_energies": "NE",
    "relativity_setting": "IREL",
    "fermi_energy_native": "EFERMI",
    "calculation_info": "INFO",
}
# Column-by-column units of the native mixed-unit records. An empty string
# marks an unknown native unit rather than assuming a physical meaning.
_SPC_COLUMN_UNITS = {
    "energy": ("degree", "eV", "arbitrary", "arbitrary", "arbitrary",
               "1", "1/angstrom", ""),
    "kx_ky": ("1/angstrom", "1/angstrom", "arbitrary", "arbitrary",
              "arbitrary", "1", ""),
}


def _exporter_version():
    try:
        return distribution_version("ase2sprkkr")
    except PackageNotFoundError:
        return None


def _spc_value(spc, field):
    value = getattr(spc, field, None)
    return value() if callable(value) else value


def _describe(dataset, long_name, units=None):
    dataset.attrs["long_name"] = long_name
    if units is not None:
        dataset.attrs["units"] = units
    return dataset


def _value(parameters, section, name):
    if parameters is None:
        return None
    group = getattr(parameters, section, None)
    option = getattr(group, name, None) if group is not None else None
    return option() if callable(option) else option


def _companion_input_file(spc):
    """Look up only names related to the recorded SPC filename."""
    source = getattr(spc, "original_filename", None) or getattr(spc, "_filename", None)
    if not source:
        return None
    source = Path(source)
    stem = source.stem
    dataset = stem.removesuffix("_ARPES_data")
    stems = {stem, dataset, dataset + "_ARPES"}
    found = {
        source.with_name(name + ext)
        for name in stems for ext in (".inp", ".in")
        if source.with_name(name + ext).is_file()
    }
    return found.pop() if len(found) == 1 else None


def _resolve_input(spc, input_parameters, input_file):
    """Explicit filename > supplied parameters > unambiguous companion file."""
    explicit = input_file is not None
    if explicit:
        candidate = Path(input_file)
        if not candidate.is_file():
            raise FileNotFoundError(f"SPR-KKR input file not found: {candidate}")
    elif input_parameters is not None:
        return input_parameters, None, input_parameters.to_string(validate=False)
    else:
        candidate = _companion_input_file(spc)

    if candidate is None:
        warnings.warn(
            "No SPR-KKR input parameters found for ARPES NeXus export; "
            "only SPC data will be exported.",
            UserWarning, stacklevel=3
        )
        return None, None, None

    from ..input_parameters.input_parameters import InputParameters
    try:
        parsed = InputParameters.from_file(candidate)
    except Exception:
        if explicit:
            raise
        warnings.warn(
            f"Cannot parse detected SPR-KKR input {candidate}; "
            "only SPC data will be exported.",
            UserWarning, stacklevel=3
        )
        return None, None, None
    return parsed, candidate, candidate.read_text(encoding="utf-8", errors="replace")


def _dataset(group, name, value):
    """Write an original numeric/string setting without inferring its meaning."""
    if value is None:
        return None
    data = np.asarray(value)
    if data.dtype.kind == "O":
        # Do not invent a representation for arbitrary Python objects.
        raise TypeError(f"Cannot serialize original input option {name}: {value!r}")
    if data.dtype.kind == "U":
        data = data.astype("S")
    return group.create_dataset(name, data=data)


def _spc_arrays(spc):
    """Use the SPC parser's named fields, shapes and native column contract."""
    mode = spc.MODE()
    columns = spc.spc_columns()
    raw = np.asarray(spc.RAW_DATA())
    if raw.ndim != 2 or raw.shape[1] != len(columns):
        raise ValueError(
            f"ARPES SPC {mode} records must have {len(columns)} columns"
        )
    fields = {
        "data": "TOTAL", "spin_up": "UP", "spin_down": "DOWN",
        "spin_polarization": "POLARIZATION", "determinant": "DETERMINANT",
    }
    if mode == "energy":
        fields.update(theta="THETA", energy="ENERGY", k_parallel="K")
    elif mode == "kx_ky":
        fields.update(k_x="KX", k_y="KY")
    else:
        raise ValueError(f"Unsupported ARPES SPC mode: {mode}")
    arrays = {
        name: np.asarray(getattr(spc, native)()) for name, native in fields.items()
    }
    shape = arrays["data"].shape
    nt, np_ = int(spc.NT()), int(spc.NP())
    if nt < 1 or np_ < 1 or any(a.shape != shape for a in arrays.values()):
        raise ValueError("Inconsistent ARPES SPC field dimensions")
    if shape[-1] != nt or (mode == "kx_ky" and shape != (np_, nt)):
        raise ValueError("ARPES SPC grid dimensions do not match NT and NP")
    if mode == "energy":
        if len(shape) == 2 and np_ != 1:
            raise ValueError("ARPES SPC NP disagrees with decoded energy-grid shape")
        if len(shape) == 3 and shape[1] != np_:
            raise ValueError("ARPES SPC NP disagrees with decoded azimuth grid")
        if len(shape) not in (2, 3):
            raise ValueError("Unexpected ARPES SPC energy-grid dimensions")
    if raw.shape != (arrays["data"].size, len(columns)):
        raise ValueError(
            f"ARPES SPC raw data must have {len(columns)} columns "
            "and match the decoded field grid"
        )
    return mode, raw, arrays


def _write_nexus(path, filename, spc, parameters, input_path, input_text, metadata):
    import h5py

    mode, raw, arrays = _spc_arrays(spc)
    shape = arrays["data"].shape
    unknown = set(metadata) - set(_EXPORT_METADATA)
    if unknown:
        raise ValueError("Unknown ARPES export metadata: " + ", ".join(sorted(unknown)))

    settings = {
        name: _value(parameters, *option)
        for name, option in _EXPORT_METADATA.items()
    }
    inpver = _value(parameters, "TASK", "INPVER")
    strver = _value(parameters, "TASK", "STRVER")
    # Old native input files can supersede the modern IP settings.
    if inpver == 0:
        settings = dict.fromkeys(settings)
    elif strver == 0:
        settings["surface_hkl"] = None
    settings.update(metadata)

    exported_at = datetime.now(timezone.utc).isoformat()
    exporter_version = _exporter_version()
    with h5py.File(path, "w") as f:
        f.attrs.update(NX_class="NXroot", default="entry", creator="ase2sprkkr",
                       file_name=filename.name, file_time=exported_at)
        if exporter_version is not None:
            f.attrs["creator_version"] = exporter_version
        entry = f.create_group("entry")
        entry.attrs.update(NX_class="NXentry", default="data")
        entry.create_dataset("title", data="SPR-KKR calculated ARPES")
        entry.create_dataset("experiment_type", data="simulation")
        program_name = entry.create_dataset("program_name", data="ase2sprkkr")
        if exporter_version is not None:
            program_name.attrs["version"] = exporter_version

        export = entry.create_group("export")
        export.attrs["NX_class"] = "NXprocess"
        export.create_dataset("program", data="ase2sprkkr")
        export.create_dataset("date", data=exported_at)
        export.create_dataset(
            "description", data="Conversion of native SPR-KKR SPC records to NeXus"
        )
        if exporter_version is not None:
            export.create_dataset("version", data=exporter_version)

        process = entry.create_group("calculation")
        process.attrs["NX_class"] = "NXprocess"
        process.create_dataset("program", data="SPR-KKR")
        process.create_dataset("description", data="One-step photoemission calculation")
        if input_text is not None:
            provenance = process.create_group("parameters")
            provenance.attrs["NX_class"] = "NXparameters"
            provenance.create_dataset("sprkkr_input", data=input_text)
            if input_path is not None:
                provenance.create_dataset("source_file", data=str(input_path))
            if parameters is not None:
                original = process.create_group("original_input")
                original.attrs["NX_class"] = "NXcollection"
                for section, names in _INPUT_OPTIONS.items():
                    options = [
                        (name, _value(parameters, section, name))
                        for name in names
                    ]
                    options = [(name, value) for name, value in options if value is not None]
                    if options:
                        group = original.create_group(section)
                        group.attrs["NX_class"] = "NXcollection"
                        for name, value in options:
                            field = _dataset(group, name, value)
                            unit = _INPUT_UNITS.get((section, name))
                            if unit is not None:
                                field.attrs["units"] = unit
                if inpver == 0 or strver == 0:
                    original.attrs["native_external_input_incomplete"] = True
        if metadata:
            supplied = process.create_group("export_metadata")
            supplied.attrs["NX_class"] = "NXcollection"
            for key, value in metadata.items():
                _dataset(supplied, key, value)

        instrument = entry.create_group("instrument")
        instrument.attrs["NX_class"] = "NXinstrument"
        source = instrument.create_group("source")
        source.attrs["NX_class"] = "NXsource"
        source.create_dataset("name", data="Virtual photon beam")
        source.create_dataset("probe", data="photon")
        if settings["photon_energy"] is not None:
            photon = _dataset(source, "energy", float(settings["photon_energy"]))
            _describe(photon, "Specified incident photon energy", "eV")
            photon.attrs["source"] = (
                "export_metadata" if "photon_energy" in metadata
                else "SPR-KKR input SPEC_PH.EPHOT"
            )
        beam = instrument.create_group("photon_geometry")
        beam.attrs["NX_class"] = "NXcollection"
        for name in ("photon_theta", "photon_phi", "photon_polarization"):
            _dataset(beam, name.removeprefix("photon_"), settings[name])
        beam.attrs["values_are"] = "specified input, not reconstructed effective parameters"
        for name, label in (
            ("theta", "Specified photon polar angle"),
            ("phi", "Specified photon azimuthal angle"),
            ("polarization", "Specified photon polarization"),
        ):
            if name in beam:
                _describe(beam[name], label,
                          "degree" if name in ("theta", "phi") else None)

        sample = entry.create_group("sample")
        sample.attrs["NX_class"] = "NXsample"
        sample.create_dataset("description", data="Calculated crystalline surface")
        if settings["surface_hkl"] is not None:
            _describe(_dataset(sample, "miller_hkl", settings["surface_hkl"]),
                      "Specified surface Miller indices", "1")
        if settings["work_function"] is not None:
            work = _dataset(sample, "work_function", settings["work_function"])
            _describe(work, "Specified work function", "eV")
            work.attrs["source"] = (
                "export_metadata" if "work_function" in metadata
                else "SPR-KKR input ENERGY.EWORK_EV"
            )

        data = entry.create_group("data")
        data.attrs.update(NX_class="NXdata", signal="data", sprkkr_mode=mode)
        data.create_dataset(
            "title", data="SPR-KKR ARPES intensity (calculated spectrum)"
        )
        main = data.create_dataset("data", data=arrays["data"])
        _describe(main, "Total ARPES intensity (arbitrary units)", "arbitrary")
        # Interpretation is a display hint, not an instruction to transpose.
        # For 3D (energy, azimuth, angle) data a default slice describes the
        # 2D image without mislabelling the 3D signal as a 2D image.
        if len(shape) == 2:
            main.attrs["interpretation"] = "image"
        elif len(shape) == 3:
            data.attrs["default_slice"] = [".", "0", "."]
        # Coordinate columns are preserved verbatim in the raw SPC below.
        if mode == "kx_ky":
            # SPC records are ordered (NP, NT): y/row index first,
            # x/column index second. Do not transpose intensity data.
            # Oblique scans have genuinely 2D coordinates; NXdata supports
            # multidimensional axis arrays via the *_indices attributes.
            # Each coordinate may span either one or both data dimensions.
            # Only collapse exactly constant coordinates: tolerances would
            # silently discard real (albeit small) deviations from a grid.
            data.attrs["axes"] = ["k_y", "k_x"]
            for name, dimension, reference, selector in (
                ("k_y", 0, arrays["k_y"][:, :1], (slice(None), 0)),
                ("k_x", 1, arrays["k_x"][:1, :], (0, slice(None))),
            ):
                coordinate = arrays[name]
                if np.allclose(coordinate, reference, rtol=0, atol=0,
                               equal_nan=True):
                    data.create_dataset(name, data=coordinate[selector])
                    data.attrs[name + "_indices"] = dimension
                else:
                    data.create_dataset(name, data=coordinate)
                    data.attrs[name + "_indices"] = [0, 1]
                _describe(
                    data[name], "Parallel photoelectron momentum " + name,
                    "1/angstrom"
                )
            if parameters is not None and _value(parameters, "SPEC_EL", "K1") is not None:
                data.attrs["scan_definition"] = "../calculation/original_input/SPEC_EL"
        elif len(shape) == 3:
            # Prefer physical axes over indices, even for nonseparable
            # coordinates. The azimuth has no physical SPC column; its
            # scan index is retained without inventing angles from input.
            energy, theta = arrays["energy"], arrays["theta"]
            data.attrs["axes"] = ["energy", "azimuth_index", "theta"]
            _describe(
                data.create_dataset("azimuth_index", data=np.arange(shape[1])),
                "Zero-based azimuth scan index (angle not present in SPC)",
                "1",
            )
            data.attrs["azimuth_index_indices"] = 1
            for name, coordinate, dim, reduced in (
                ("energy", energy, 0, energy[:, :1, :1]),
                ("theta", theta, 2, theta[:1, :1, :]),
            ):
                if np.allclose(coordinate, reduced, rtol=0, atol=0, equal_nan=True):
                    selector = (slice(None), 0, 0) if dim == 0 else (0, 0, slice(None))
                    data.create_dataset(name, data=coordinate[selector])
                    data.attrs[name + "_indices"] = dim
                else:
                    data.create_dataset(name, data=coordinate)
                    data.attrs[name + "_indices"] = [0, 1, 2]
        else:
            energy, theta = arrays["energy"], arrays["theta"]
            data.attrs["axes"] = ["energy", "theta"]
            for name, coordinate, dim, reduced, selector in (
                ("energy", energy, 0, energy[:, :1], (slice(None), 0)),
                ("theta", theta, 1, theta[:1, :], (0, slice(None))),
            ):
                if np.allclose(coordinate, reduced, rtol=0, atol=0, equal_nan=True):
                    data.create_dataset(name, data=coordinate[selector])
                    data.attrs[name + "_indices"] = dim
                else:
                    data.create_dataset(name, data=coordinate)
                    data.attrs[name + "_indices"] = [0, 1]

        if mode == "energy":
            _describe(data["energy"], "Initial-state energy E - E_F", "eV")
            _describe(data["theta"], "Photoelectron polar angle", "degree")
            _describe(
                data.create_dataset("k_parallel", data=arrays["k_parallel"]),
                "Parallel photoelectron momentum", "1/angstrom",
            )
            data.attrs["k_parallel_indices"] = list(range(len(shape)))
        for name in ("spin_up", "spin_down", "spin_polarization", "determinant"):
            data.create_dataset(name, data=arrays[name])
        _describe(data["spin_up"], "Spin-up ARPES intensity (arbitrary units)",
                  "arbitrary")
        _describe(data["spin_down"], "Spin-down ARPES intensity (arbitrary units)",
                  "arbitrary")
        _describe(data["spin_polarization"], "Spin polarization", "1")
        data["spin_polarization"].attrs["formula"] = "(I_1 - I_2) / (I_1 + I_2)"
        _describe(data["determinant"], "Native SPR-KKR spectral determinant")
        # The determinant's physical units are not documented by the SPC
        # parser, so do not assign a potentially incorrect unit.
        data.attrs["auxiliary_signals"] = [
            "spin_up", "spin_down", "spin_polarization", "determinant"
        ]

        # Add HDF5 dimension labels/scales as a compatibility aid for generic
        # HDF5 viewers. NeXus axes and *_indices remain authoritative. A
        # coupled (multidimensional) coordinate is never flattened to invent
        # a one-dimensional scale.
        for dimension, axis_name in enumerate(data.attrs["axes"]):
            main.dims[dimension].label = axis_name
            coordinate = data[axis_name]
            if coordinate.ndim == 1 and coordinate.shape[0] == shape[dimension]:
                coordinate.make_scale(axis_name)
                main.dims[dimension].attach_scale(coordinate)

        native = entry.create_group("sprkkr_spc")
        native.attrs["NX_class"] = "NXcollection"
        native.create_dataset("mode", data=mode)
        spc_name = (
            getattr(spc, "original_filename", None)
            or getattr(spc, "_filename", None)
        )
        if spc_name is not None:
            native.create_dataset("source_filename", data=Path(spc_name).name)
        raw_ds = native.create_dataset("raw_data", data=raw)
        raw_ds.attrs["columns"] = spc.spc_columns()
        raw_ds.attrs["column_units"] = _SPC_COLUMN_UNITS[mode]
        raw_ds.attrs["long_name"] = "Unmodified numerical records from SPR-KKR SPC"
        raw_ds.attrs["record_shape"] = shape
        raw_ds.attrs["record_order"] = (
            "NE, NP, NT (NT fastest)" if len(shape) == 3 else
            ("NE, NT (NT fastest)" if mode == "energy"
             else "NP, NT (NT fastest)")
        )
        native.create_dataset("nt", data=int(spc.NT()))
        native.create_dataset("np", data=int(spc.NP()))
        native.create_dataset("comment", data=str(spc.COMMENT()))
        header = [
            (name, _spc_value(spc, field))
            for name, field in _SPC_HEADER.items()
        ]
        header = [(name, value) for name, value in header if value is not None]
        if header:
            original_header = native.create_group("header")
            original_header.attrs["NX_class"] = "NXcollection"
            for name, value in header:
                ds = _dataset(original_header, name, value)
                if name in ("number_of_sites", "number_of_atom_types",
                            "number_of_energies"):
                    ds.attrs["units"] = "1"
                if name == "fermi_energy_native":
                    ds.attrs["long_name"] = (
                        "Native SPC EFERMI value (unit not inferred)"
                    )


def write_arpes_nexus(spc, filename, *, input_parameters=None,
                      input_file=None, metadata=None):
    """Export native SPC without inferring unavailable instrument metadata.

    An explicit input_file takes precedence; missing input issues UserWarning.
    Write to a temporary path and replace the target only on success.
    """
    try:
        import h5py  # noqa: F401
    except ImportError as exc:
        raise ImportError("ARPES NeXus export needs h5py (ase2sprkkr[nexus])") from exc

    parameters, source_path, input_text = _resolve_input(
        spc, input_parameters, input_file
    )
    target = Path(filename)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            prefix="." + target.name + ".", suffix=".tmp",
            dir=target.parent, delete=False
        ) as tmp:
            temporary = Path(tmp.name)
        _write_nexus(temporary, target, spc, parameters, source_path,
                     input_text, dict(metadata or {}))
        os.replace(temporary, target)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return target
