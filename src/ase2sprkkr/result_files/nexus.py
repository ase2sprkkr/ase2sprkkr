"""Export native SPR-KKR ARPES SPC data as NeXus/HDF5.

The SPC file owns sampled intensities and coordinates. InputParameters owns
the *specified* experiment settings. This module does not recompute KKRSPEC
spin, photon, or reciprocal-space transformations.
"""

from datetime import datetime, timezone
from pathlib import Path
import os
import tempfile
import warnings

import numpy as np


# A single declarative map for optional, *original* SPR-KKR settings.
# Values in this map are not inferred effective values.
_INPUT_OPTIONS = {
    "TASK": ("TASK", "INPVER", "STRVER", "MILLER_HKL", "CRYS_VEC"),
    "ENERGY": ("EMIN", "EMAX", "NE", "EWORK_EV", "IMV_INI_EV", "IMV_FIN_EV"),
    "SPEC_PH": ("EPHOT", "THETA", "PHI", "POL_P", "NPOL", "ICIRC",
                "IDREH", "ALQ", "DELQ", "IFSP", "THETA_FIX"),
    "SPEC_EL": ("THETA", "PHI", "NT", "NP", "KA", "K1", "K2",
                "NK1", "NK2", "SPOL", "PSPIN", "POL_E",
                "Q1", "Q2", "Q3", "Q4", "BETA1", "BETA2", "ROTAXIS"),
    "SPEC_STR": ("SURF_BAR", "N_LAYDBL", "N_LAYER"),
}
_EXPORT_METADATA = {
    "photon_energy": ("SPEC_PH", "EPHOT"),
    "photon_theta": ("SPEC_PH", "THETA"),
    "photon_phi": ("SPEC_PH", "PHI"),
    "photon_polarization": ("SPEC_PH", "POL_P"),
    "work_function": ("ENERGY", "EWORK_EV"),
    "surface_hkl": ("TASK", "MILLER_HKL"),
}


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
    raw = np.asarray(spc.RAW_DATA())
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

    with h5py.File(path, "w") as f:
        f.attrs.update(NX_class="NXroot", default="entry", creator="ase2sprkkr",
                       file_name=filename.name,
                       file_time=datetime.now(timezone.utc).isoformat())
        entry = f.create_group("entry")
        entry.attrs.update(NX_class="NXentry", default="data")
        entry.create_dataset("title", data="SPR-KKR calculated ARPES")
        entry.create_dataset("experiment_type", data="simulation")

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
                        for name, value in options:
                            _dataset(group, name, value)
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
            photon.attrs["units"] = "eV"
        beam = instrument.create_group("photon_geometry")
        beam.attrs["NX_class"] = "NXcollection"
        for name in ("photon_theta", "photon_phi", "photon_polarization"):
            _dataset(beam, name.removeprefix("photon_"), settings[name])
        beam.attrs["values_are"] = "specified input, not reconstructed effective parameters"
        for name in ("theta", "phi"):
            if name in beam:
                beam[name].attrs["units"] = "degree"

        sample = entry.create_group("sample")
        sample.attrs["NX_class"] = "NXsample"
        sample.create_dataset("description", data="Calculated crystalline surface")
        if settings["surface_hkl"] is not None:
            _dataset(sample, "miller_hkl", settings["surface_hkl"])

        data = entry.create_group("data")
        data.attrs.update(NX_class="NXdata", signal="data", sprkkr_mode=mode)
        main = data.create_dataset("data", data=arrays["data"])
        main.attrs["units"] = "arbitrary_units"
        # Coordinate columns are preserved verbatim in the raw SPC below.
        if mode == "kx_ky":
            # NeXus regular axes cannot encode an oblique reciprocal-space
            # basis. Display its projections on orthogonal kx/ky axes; native
            # 2D coordinates are retained exactly in sprkkr_spc/raw_data.
            warning = (
                "ARPES kx_ky coordinates are represented as orthogonal "
                "k_x/k_y axes in NXdata; oblique scans may be displayed "
                "with distorted geometry. Exact coordinates are in "
                "entry/sprkkr_spc/raw_data (columns KX, KY)."
            )
            warnings.warn(warning, UserWarning, stacklevel=3)
            data.attrs["coordinate_warning"] = warning
            data.attrs["axes"] = ["k_y", "k_x"]
            data.attrs["k_y_indices"] = 0
            data.attrs["k_x_indices"] = 1
            kx = data.create_dataset("k_x", data=arrays["k_x"][0, :])
            ky = data.create_dataset("k_y", data=arrays["k_y"][:, 0])
            kx.attrs["units"] = ky.attrs["units"] = "1/angstrom"
            if parameters is not None and _value(parameters, "SPEC_EL", "K1") is not None:
                data.attrs["scan_definition"] = "../calculation/original_input/SPEC_EL"
        elif len(shape) == 3:
            data.attrs["axes"] = ["energy_index", "azimuth_index", "angle_index"]
            for axis, count, index in zip(
                ("energy_index", "azimuth_index", "angle_index"), shape, range(3)
            ):
                data.create_dataset(axis, data=np.arange(count))
                data.attrs[axis + "_indices"] = index
            for coord in ("theta", "energy"):
                data.create_dataset(coord, data=arrays[coord])
                data.attrs[coord + "_indices"] = [0, 1, 2]
        else:
            energy, theta = arrays["energy"], arrays["theta"]
            separable = (
                np.allclose(energy, energy[:, :1], equal_nan=True)
                and np.allclose(theta, theta[:1, :], equal_nan=True)
            )
            if separable:
                data.attrs["axes"] = ["energy", "theta"]
                data.attrs["energy_indices"] = 0
                data.attrs["theta_indices"] = 1
                data.create_dataset("energy", data=energy[:, 0])
                data.create_dataset("theta", data=theta[0, :])
            else:
                data.attrs["axes"] = ["energy_index", "angle_index"]
                for axis, count, index in zip(
                    ("energy_index", "angle_index"), shape, range(2)
                ):
                    data.create_dataset(axis, data=np.arange(count))
                    data.attrs[axis + "_indices"] = index
                for coord in ("energy", "theta"):
                    data.create_dataset(coord, data=arrays[coord])
                    data.attrs[coord + "_indices"] = [0, 1]

        if mode == "energy":
            data["energy"].attrs["units"] = "eV"
            data["theta"].attrs["units"] = "degree"
            data.create_dataset("k_parallel", data=arrays["k_parallel"])
            data["k_parallel"].attrs["units"] = "1/angstrom"
            data.attrs["k_parallel_indices"] = list(range(len(shape)))
        for name in ("spin_up", "spin_down", "spin_polarization", "determinant"):
            data.create_dataset(name, data=arrays[name])
        data["spin_polarization"].attrs["units"] = "1"
        data["spin_polarization"].attrs["formula"] = "(I_1 - I_2) / (I_1 + I_2)"
        data.attrs["auxiliary_signals"] = [
            "spin_up", "spin_down", "spin_polarization", "determinant"
        ]

        native = entry.create_group("sprkkr_spc")
        native.attrs["NX_class"] = "NXcollection"
        raw_ds = native.create_dataset("raw_data", data=raw)
        raw_ds.attrs["columns"] = spc.spc_columns()
        native.create_dataset("nt", data=int(spc.NT()))
        native.create_dataset("np", data=int(spc.NP()))
        native.create_dataset("comment", data=str(spc.COMMENT()))


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
