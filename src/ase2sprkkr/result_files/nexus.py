"""NeXus export of simulated SPR-KKR photoemission spectra.

This writes NeXus base classes. It deliberately does not claim compliance
with NXmpes_arpes, which requires a physical energy-dispersion scheme.
"""

from datetime import datetime, timezone
from pathlib import Path
import warnings

import numpy as np


# Native SPR-KKR input-to-NeXus provenance mapping. Explicit metadata overrides
# are allowed for inputs not preserved with an independently loaded SPC file.
_NATIVE_PARAMETERS = {
    "photon_energy": ("SPEC_PH", "EPHOT"),
    "work_function": ("ENERGY", "EWORK_EV"),
    "photon_theta": ("SPEC_PH", "THETA"),
    "photon_phi": ("SPEC_PH", "PHI"),
    "photon_polarization": ("SPEC_PH", "POL_P"),
    "spin_mode": ("SPEC_EL", "SPOL"),
    "spin_axis": ("SPEC_EL", "PSPIN"),
    "surface_hkl": ("TASK", "MILLER_HKL"),
    "scan_k1": ("SPEC_EL", "K1"),
    "scan_k2": ("SPEC_EL", "K2"),
    "scan_ka": ("SPEC_EL", "KA"),
    "beta1": ("SPEC_EL", "BETA1"),
    "beta2": ("SPEC_EL", "BETA2"),
    "rotaxis": ("SPEC_EL", "ROTAXIS"),
    "electron_theta": ("SPEC_EL", "THETA"),
    "electron_phi": ("SPEC_EL", "PHI"),
    "scan_nt": ("SPEC_EL", "NT"),
    "scan_np": ("SPEC_EL", "NP"),
    "electron_typ": ("SPEC_EL", "TYP"),
    "photon_alq": ("SPEC_PH", "ALQ"),
    "photon_delq": ("SPEC_PH", "DELQ"),
    "photon_npol": ("SPEC_PH", "NPOL"),
    "photon_icirc": ("SPEC_PH", "ICIRC"),
    "photon_idreh": ("SPEC_PH", "IDREH"),
    "surface_barrier": ("SPEC_STR", "SURF_BAR"),
    "initial_damping": ("ENERGY", "IMV_INI_EV"),
    "final_damping": ("ENERGY", "IMV_FIN_EV"),
    "photon_ifsp": ("SPEC_PH", "IFSP"),
    "photon_theta_fix": ("SPEC_PH", "THETA_FIX"),
    "spin_polarization": ("SPEC_EL", "POL_E"),
    "spin_pol0": ("SPEC_EL", "POL0"),
    "spin_pol0l": ("SPEC_EL", "POL0L"),
    "spin_q1": ("SPEC_EL", "Q1"),
    "spin_q2": ("SPEC_EL", "Q2"),
    "spin_q3": ("SPEC_EL", "Q3"),
    "spin_q4": ("SPEC_EL", "Q4"),
    "input_version": ("TASK", "INPVER"),
    "structure_version": ("TASK", "STRVER"),
    "external_input_filename": ("TASK", "INPFIL"),
    "external_structure_filename": ("TASK", "STRFIL"),
}


def _photon_setup(values):
    """Effective parameters after SPEC_INPUT polarization correction.

    See SPR-KKR 9.7 spec_inputs.f (POL_P dispatch and ICIRC correction).
    POL_P overrides the low-level input values, and IFSP is read only
    when POL_P was supplied.
    """
    values = dict(values)
    polar = values["photon_polarization"]
    npol = 1 if values["photon_npol"] is None else int(values["photon_npol"])
    icirc = 1 if values["photon_icirc"] is None else int(values["photon_icirc"])
    idreh = 0 if values["photon_idreh"] is None else int(values["photon_idreh"])
    alq = 45. if values["photon_alq"] is None else float(values["photon_alq"])
    delq = 0. if values["photon_delq"] is None else float(values["photon_delq"])
    if polar is not None:
        presets = {"P": (1, 1, 1), "S": (2, 1, 1),
                   "C+": (1, 2, 1), "C-": (2, 2, 1)}
        key = str(polar).strip().upper()
        if key not in presets:
            raise ValueError("Unknown SPEC_PH.POL_P polarization: " + str(polar))
        npol, icirc, idreh = presets[key]
    else:
        # Without POL_P, the native SPEC_INPUT does not even read IFSP.
        values["photon_ifsp"] = 0
        values["photon_theta_fix"] = None
    if npol == 0:
        icirc = 1
    if icirc in (1, 4):
        delq, idreh = 0., 0
        if npol in (0, 1):
            alq = 0.
        elif npol == 2:
            alq = 90.
    elif icirc in (2, 3):
        alq, delq = 45., 90.
        idreh = 1 if npol in (1, 3) else -1
    elif idreh:
        idreh = 1 if idreh > 0 else -1
    values.update(photon_npol=npol, photon_icirc=icirc,
                  photon_idreh=idreh, photon_alq=alq, photon_delq=delq)
    if values["photon_ifsp"] == 1 and values["photon_theta_fix"] is None:
        values["photon_theta_fix"] = values["photon_theta"]
    return values


def _spin_amplitudes(values, *, use_polarization):
    """Native Q1..Q4 choices from SPEC_INPUT; not a detector efficiency."""
    amplitudes = [1.+0j, 0j, 0j, 1.+0j]
    polarization = values["spin_polarization"]
    presets = {
        "PZ": (1, 0, 0, 1),
        "PX": (.707, .707, .707, -.707),
        "PY": (.707, .707j, .707, -.707j),
    }
    if use_polarization and polarization is not None:
        amplitudes = presets.get(str(polarization).upper(), amplitudes)
    else:
        amplitudes = tuple(
            amplitudes[i] if values["spin_q" + str(i+1)] is None
            else complex(values["spin_q" + str(i+1)])
            for i in range(4)
        )
    return np.asarray(amplitudes, dtype=complex)


def _companion_input_file(spc):
    """Find an unambiguous input next to the SPC source; never guess among ties."""
    source = (
        getattr(spc, "original_filename", None)
        or getattr(spc, "_source_filename", None)
    )
    if source is None:
        return None
    source = Path(source)
    if not source.name:
        return None
    stem = source.stem
    dataset = stem[:-len("_ARPES_data")] if stem.endswith("_ARPES_data") else stem
    stems = {stem, dataset, dataset + "_ARPES"}
    directory = source.parent
    matching = {
        directory / (name + ext)
        for name in stems for ext in (".inp", ".in")
        if (directory / (name + ext)).is_file()
    }
    if len(matching) == 1:
        return matching.pop()
    if matching:
        return None
    nearby = {
        path for ext in (".inp", ".in")
        for path in directory.glob("*" + ext) if path.is_file()
    }
    return next(iter(nearby)) if len(nearby) == 1 else None


def _input_for_export(spc, input_parameters, input_file):
    """Resolve input with explicit file > given object > companion discovery."""
    if input_file is not None:
        source = Path(input_file)
        if not source.is_file():
            raise FileNotFoundError(f"SPR-KKR input file not found: {source}")
    elif input_parameters is not None:
        return input_parameters, None
    else:
        source = _companion_input_file(spc)

    if source is None:
        warnings.warn(
            "No SPR-KKR input parameters found for ARPES NeXus export; "
            "the spectrum is exported without a verified calculation setup. "
            "Pass input_file=... or input_parameters=... to include it.",
            UserWarning, stacklevel=3
        )
        return None, None

    from ..input_parameters.input_parameters import InputParameters
    try:
        return InputParameters.from_file(source), source
    except Exception:
        # A deliberately specified input must fail rather than silently
        # being replaced with incorrect calculated setup.
        if input_file is not None:
            raise
        warnings.warn(
            f"Could not read detected SPR-KKR input {source}; "
            "exporting spectrum without calculation setup.",
            UserWarning, stacklevel=3
        )
        return None, None

def write_arpes_nexus(spc, filename, *, input_parameters=None,
                      input_file=None, metadata=None):
    """Write an ARPES SPC result to NeXus/HDF5 and return the output path.

    The original eight data columns are preserved. The energy origin is the
    Fermi level; k_parallel is in inverse angstrom, as verified in
    SPR-KKR 9.7 spec_sumup_int.f. This is an energy/angle ARPES scan.
    """
    try:
        import h5py
    except ImportError as exc:
        raise ImportError("ARPES NeXus export requires h5py (install ase2sprkkr[nexus])") from exc

    mode = spc.MODE() if hasattr(spc, "MODE") else "energy"
    if mode not in ("energy", "kx_ky"):
        raise ValueError(f"Unknown ARPES SPC mode: {mode}")
    values = {
        "total": spc.TOTAL(),
        "spin_up": spc.UP(),
        "spin_down": spc.DOWN(),
        "spin_polarization": spc.POLARIZATION(),
        "determinant": spc.DETERMINANT(),
    }
    if mode == "energy":
        values.update(theta=spc.THETA(), energy=spc.ENERGY(), k_parallel=spc.K())
    else:
        values.update(k_x=spc.KX(), k_y=spc.KY())
    arrays = {name: np.asarray(value) for name, value in values.items()}
    if mode == "energy" and int(spc.NP()) != 1:
        raise ValueError("Energy/angle SPC requires NP=1")
    shape = arrays["total"].shape
    if len(shape) != 2 or min(shape) == 0:
        raise ValueError("ARPES spectrum must form a non-empty (energy, angle) grid")
    if shape[1] != int(spc.NT()) or (mode == "kx_ky" and shape[0] != int(spc.NP())):
        raise ValueError("ARPES SPC dimensions do not match the declared grid")
    if any(value.shape != shape for value in arrays.values()):
        raise ValueError("ARPES SPC columns have inconsistent dimensions")
    if any(not np.issubdtype(value.dtype, np.number) for value in arrays.values()):
        raise ValueError("ARPES SPC arrays must be numeric")

    if mode == "energy":
        energy, theta = arrays["energy"], arrays["theta"]
        rectilinear = (
            np.allclose(energy, energy[:, :1], equal_nan=True)
            and np.allclose(theta, theta[:1, :], equal_nan=True)
        )
    else:
        kx, ky = arrays["k_x"], arrays["k_y"]
        # SPR-KKR generates a regular reciprocal-space mesh with JP and JT
        # as its two grid indices. Use the full spans to reduce the effect
        # of limited precision in the text SPC coordinates.
        k_origin = np.array([kx[0, 0], ky[0, 0]], dtype=float)
        k_basis = np.zeros((2, 2), dtype=float)
        if shape[0] > 1:
            k_basis[0] = (
                np.array([kx[-1, 0], ky[-1, 0]]) - k_origin
            ) / (shape[0] - 1)
        if shape[1] > 1:
            k_basis[1] = (
                np.array([kx[0, -1], ky[0, -1]]) - k_origin
            ) / (shape[1] - 1)

    raw_data = np.asarray(spc.RAW_DATA())
    ncols = 8 if mode == "energy" else 7
    if raw_data.ndim != 2 or raw_data.shape != (shape[0] * shape[1], ncols):
        raise ValueError(f"ARPES SPC {mode} raw data must have {ncols} columns and match the grid")
    input_parameters, resolved_input_file = _input_for_export(
        spc, input_parameters, input_file
    )
    input_text = (
        input_parameters.to_string(validate=False)
        if input_parameters is not None else None
    )
    metadata = dict(metadata or {})
    unknown = set(metadata) - set(_NATIVE_PARAMETERS)
    if unknown:
        raise ValueError("Unknown ARPES export metadata: " + ", ".join(sorted(unknown)))

    def value_from_input(section, name):
        if input_parameters is None:
            return None
        group = getattr(input_parameters, section, None)
        option = getattr(group, name, None) if group is not None else None
        return option() if callable(option) else option

    derived = {
        key: value_from_input(*source)
        for key, source in _NATIVE_PARAMETERS.items()
    }
    input_version = derived["input_version"]
    structure_version = derived["structure_version"]
    if input_version == 0:
        # The old input.inp/barrier file determines the actual photoemission
        # setup, not the values retained in modern InputParameters.
        filenames = {key: derived[key] for key in (
            "external_input_filename", "external_structure_filename"
        )}
        derived = {key: None for key in derived}
        derived.update(input_version=0, structure_version=structure_version)
        derived.update(filenames)
    elif structure_version == 0:
        # Do not claim that the Python Miller indices are the actual surface.
        derived["surface_hkl"] = None
    effective = {**derived, **metadata}
    # Resolve the polarization as KKRSPEC does, after applying overrides.
    # Standalone SPC exports without photon metadata must not acquire
    # fabricated instrument settings.
    polarization_keys = {
        "photon_polarization", "photon_npol", "photon_icirc",
        "photon_idreh", "photon_alq", "photon_delq", "photon_ifsp"
    }
    if (input_version != 0 and input_parameters is not None) or (
        polarization_keys.intersection(metadata)
    ):
        effective = _photon_setup(effective)
    photon_energy = effective["photon_energy"]
    work_function = effective["work_function"]
    photon_theta = effective["photon_theta"]
    photon_phi = effective["photon_phi"]
    photon_polarization = effective["photon_polarization"]
    spin_mode = effective["spin_mode"]
    spin_axis = effective["spin_axis"] if spin_mode == 4 else None
    surface_hkl = effective["surface_hkl"]
    scan_k1 = effective["scan_k1"]
    scan_k2 = effective["scan_k2"]
    scan_ka = effective["scan_ka"]
    beta1 = effective["beta1"]
    beta2 = effective["beta2"]
    rotaxis = effective["rotaxis"]
    photon_ifsp = effective["photon_ifsp"]
    photon_theta_fix = effective["photon_theta_fix"]
    spin_polarization = effective["spin_polarization"]
    # For ordinary ARPES, POL0/POL0L are only read by SPEC_INPUT on the
    # SPLEED branch (NSPLEE=0), not as active spin analyzer axes.
    spin_q_input = {
        key: effective[key] for key in ("spin_q1", "spin_q2", "spin_q3", "spin_q4")
    }
    if input_text is not None:
        import re
        polarization_in_input = bool(re.search(
            r"(?mi)^\s*POL_E\s*(?:=|\s)\s*\w+", input_text
        ))
    else:
        polarization_in_input = False
    use_spin_polarization = (
        input_version != 0 and
        (polarization_in_input or "spin_polarization" in metadata)
    )
    spin_amplitudes = _spin_amplitudes(effective, use_polarization=use_spin_polarization)

    electron_theta = effective["electron_theta"]
    electron_phi = effective["electron_phi"]
    scan_nt = effective["scan_nt"]
    scan_np = effective["scan_np"]
    electron_typ = effective["electron_typ"]
    photon_alq = effective["photon_alq"]
    photon_delq = effective["photon_delq"]
    photon_npol = effective["photon_npol"]
    photon_icirc = effective["photon_icirc"]
    photon_idreh = effective["photon_idreh"]
    surface_barrier = effective["surface_barrier"]
    initial_damping = effective["initial_damping"]
    final_damping = effective["final_damping"]
    if spin_mode == 4 and spin_axis is None and input_parameters is not None:
        # SPEC_INPUT defaults to +z if no explicit PSPIN is supplied.
        spin_axis = (0., 0., 1.)

    filename = Path(filename)
    with h5py.File(filename, "w") as file:
        file.attrs["NX_class"] = "NXroot"
        file.attrs["default"] = "entry"
        file.attrs["file_name"] = filename.name
        file.attrs["h5py_version"] = h5py.__version__
        file.attrs["file_time"] = datetime.now(timezone.utc).isoformat()
        file.attrs["creator"] = "ase2sprkkr"
        entry = file.create_group("entry")
        entry.attrs.update(NX_class="NXentry", default="data")
        entry.create_dataset("title", data="SPR-KKR calculated ARPES spectrum")
        entry.create_dataset("experiment_type", data="simulation")

        instrument = entry.create_group("instrument")
        instrument.attrs["NX_class"] = "NXinstrument"
        instrument.create_dataset("name", data="SPR-KKR virtual ARPES spectrometer")
        analyzer = instrument.create_group("electronanalyzer")
        analyzer.attrs["NX_class"] = "NXelectronanalyzer"
        analyzer.create_dataset(
            "description",
            data=("Ideal numerical angle/energy sampling of the calculated photocurrent; "
                  "no physical electron analyzer, energy dispersion, angular acceptance "
                  "or instrumental energy/angular resolution is modeled."),
        )
        photon_source = instrument.create_group("source")
        photon_source.attrs["NX_class"] = "NXsource"
        photon_source.create_dataset("name", data="Virtual incident photon beam")
        photon_source.create_dataset("probe", data="photon")
        if photon_energy is not None:
            d = photon_source.create_dataset("energy", data=float(photon_energy))
            d.attrs["units"] = "eV"
        photon_geometry = instrument.create_group("photon_geometry")
        photon_geometry.attrs["NX_class"] = "NXcollection"
        photon_geometry.create_dataset("description", data=(
            "SPR-KKR SPEC_PH polar/azimuth convention; theta=180 degrees "
            "is normal incidence. No laboratory-frame transformation is asserted."
        ))
        photon_geometry.create_dataset("ifsp", data=int(photon_ifsp or 0))
        if photon_ifsp == 1:
            if photon_theta_fix is None:
                raise ValueError("IFSP=1 requires THETA_FIX or the photon THETA")
            fixed = photon_geometry.create_dataset("theta_fix", data=float(photon_theta_fix))
            fixed.attrs["units"] = "degree"
            photon_geometry.create_dataset("geometry_mode", data="emission_angle_dependent")
            # In SPEC_VBRUN, THEQ = THETA_electron - FIXTHQ, then SPEC_TILT.
            # This is the pre-tilt angle, not a fabricated lab-frame angle.
            if mode == "energy":
                angle = photon_geometry.create_dataset(
                    "polar_angle_before_tilt",
                    data=arrays["theta"] - float(photon_theta_fix)
                )
                angle.attrs["units"] = "degree"
                angle.attrs["long_name"] = "Native photon polar angle before SPEC_TILT"
            else:
                photon_geometry.create_dataset(
                    "angular_grid_note",
                    data="The constant-energy SPC does not store emission-angle grid."
                )
            photon_geometry.create_dataset("calculation", data=(
                "SPEC_VBRUN: THEQ = THETA_electron - FIXTHQ; "
                "SPEC_TILT may rotate photon incidence for sample tilt"
            ))
        else:
            photon_geometry.create_dataset("geometry_mode", data="fixed_incidence")
            if photon_theta is not None:
                d = photon_geometry.create_dataset("theta", data=float(photon_theta))
                d.attrs["units"] = "degree"
        if photon_phi is not None:
            d = photon_geometry.create_dataset("phi", data=float(photon_phi))
            d.attrs["units"] = "degree"
        if photon_polarization is not None:
            photon_geometry.create_dataset("polarization", data=str(photon_polarization))
        for name, value, units in (
            ("alq", photon_alq, "degree"), ("delq", photon_delq, "degree"),
            ("npol", photon_npol, None), ("icirc", photon_icirc, None),
            ("idreh", photon_idreh, None)
        ):
            if value is not None:
                d = photon_geometry.create_dataset(name, data=value)
                if units is not None:
                    d.attrs["units"] = units

        sample = entry.create_group("sample")
        sample.attrs["NX_class"] = "NXsample"
        sample.create_dataset("description", data="Computational crystalline surface")
        if surface_hkl is not None:
            orient = sample.create_group("surface_orientation")
            orient.attrs["NX_class"] = "NXcollection"
            orient.create_dataset("miller_hkl", data=np.asarray(surface_hkl, dtype=int))
            orient.create_dataset("description", data=(
                "Miller indices provided to SPR-KKR; do not imply an absolute "
                "laboratory orientation or validated coordinate transformation."
            ))

        electron_sampling = analyzer.create_group("sampling")
        electron_sampling.attrs["NX_class"] = "NXcollection"
        electron_sampling.create_dataset("description", data=(
            "SPR-KKR native electron angular/k scan. Emission angles and "
            "k_parallel are stored exactly as output by KKRSPEC."
        ))
        if scan_k1 is not None:
            electron_sampling.create_dataset("mode", data="parallel_momentum")
            electron_sampling.create_dataset("k1", data=np.asarray(scan_k1, dtype=float))
            electron_sampling["k1"].attrs["units"] = "1/angstrom"
            # SPEC_INPUT produces its one-dimensional K mesh from KA and K1.
            # K1 is a displacement, not the end-point coordinate.
            start = np.zeros(2) if scan_ka is None else np.asarray(scan_ka, dtype=float)
            step = np.asarray(scan_k1, dtype=float)
            if start.shape == (2,) and step.shape == (2,):
                length = np.linalg.norm(step)
                if shape[1] == 1:
                    kmesh = start.reshape(1, 2)
                else:
                    fraction = np.arange(shape[1], dtype=float) / (shape[1] - 1)
                    # Match the actual SPEC_INPUT formula, including its
                    # small denominator regularization.
                    displacement = length * step / np.sqrt(length**2 + 1e-10)
                    kmesh = start[None, :] + fraction[:, None] * displacement[None, :]
                for i, name in enumerate(("k_x", "k_y")):
                    ds = electron_sampling.create_dataset(name, data=kmesh[:, i])
                    ds.attrs["units"] = "1/angstrom"
                    ds.attrs["long_name"] = "Native 1D reciprocal-space sampling coordinate"
                electron_sampling.create_dataset(
                    "coordinate_reference", data="k_x/k_y include the KA offset"
                )
            for name, value in (("k2", scan_k2), ("ka", scan_ka)):
                if value is not None:
                    ds = electron_sampling.create_dataset(name, data=np.asarray(value, dtype=float))
                    ds.attrs["units"] = "1/angstrom"
        else:
            electron_sampling.create_dataset("mode", data=(
                "angular" if input_parameters is not None else "unspecified"
            ))
        if beta1 is not None or beta2 is not None:
            tilt = sample.create_group("tilt")
            tilt.attrs["NX_class"] = "NXcollection"
            for name, value in (("beta1", beta1), ("beta2", beta2)):
                if value is not None:
                    ds = tilt.create_dataset(name, data=float(value))
                    ds.attrs["units"] = "degree"
            if rotaxis is not None:
                tilt.create_dataset("rotation_axis", data=np.asarray(rotaxis, dtype=float))
        for name, value in (
            ("theta_range", electron_theta), ("phi_range", electron_phi),
            ("nt", scan_nt), ("np", scan_np), ("typ", electron_typ)
        ):
            if value is not None:
                d = electron_sampling.create_dataset(name, data=np.asarray(value))
                if name.endswith("_range"):
                    d.attrs["units"] = "degree"
        if surface_barrier is not None:
            sample.create_dataset("barrier_parameters", data=np.asarray(surface_barrier, dtype=float))
        if spin_mode is not None:
            projection = instrument.create_group("spin_projection")
            projection.attrs["NX_class"] = "NXcollection"
            projection.create_dataset("spol", data=int(spin_mode))
            projection.create_dataset("description", data=(
                "SPOL=1 unpolarized, 2 spin channels, 4 spin density matrix. "
                "For SPOL=4, the native (channel1-channel2) sign is opposite "
                "to the dot product of spin polarization with PSPIN."
            ))
            if spin_axis is not None:
                projection.create_dataset("pspin", data=np.asarray(spin_axis, dtype=float))
            projection.create_dataset("q_amplitudes", data=spin_amplitudes)
            projection["q_amplitudes"].attrs["long_name"] = (
                "Effective native Q1..Q4 complex spin amplitudes"
            )
            if use_spin_polarization and spin_polarization is not None:
                projection.create_dataset("pol_e", data=str(spin_polarization))
            projection.create_dataset("coordinate_convention", data=(
                "POL0 and POL0L spin-axis vectors are consumed by SPR-KKR "
                "on the SPLEED branch, not for ARPES (NSPLEE != 0)."
            ))

        process = entry.create_group("calculation")
        process.attrs["NX_class"] = "NXprocess"
        process.create_dataset("program", data="SPR-KKR")
        process.create_dataset("description", data="Theoretical one-step photoemission (kkrspec)")
        if input_version == 0 or structure_version == 0:
            external = process.create_group("external_native_input")
            external.attrs["NX_class"] = "NXcollection"
            external.create_dataset("geometry_complete", data=False)
            external.create_dataset("description", data=(
                "The run uses external native SPR-KKR input and/or structure "
                "files. Their physical metadata cannot be inferred solely "
                "from InputParameters; unset fields are not asserted."
            ))
            external.create_dataset("input_version", data=(
                -1 if input_version is None else int(input_version)
            ))
            external.create_dataset("structure_version", data=(
                -1 if structure_version is None else int(structure_version)
            ))
            for name, value in (
                ("input_filename", derived.get("external_input_filename")),
                ("structure_filename", derived.get("external_structure_filename"))
            ):
                if value is not None:
                    external.create_dataset(name, data=str(value))
        if metadata:
            supplied = process.create_group("export_metadata")
            supplied.attrs["NX_class"] = "NXcollection"
            supplied.create_dataset("keys_overridden_or_supplied", data=sorted(metadata))
            for name, value in sorted(metadata.items()):
                if value is not None:
                    supplied.create_dataset(name, data=np.asarray(value) if not isinstance(value, str) else value)
        if any(value is not None for value in (initial_damping, final_damping)):
            damping = process.create_group("imaginary_potential")
            damping.attrs["NX_class"] = "NXcollection"
            for name, value in (("initial_state", initial_damping), ("final_state", final_damping)):
                if value is not None:
                    d = damping.create_dataset(name, data=float(value))
                    d.attrs["units"] = "eV"
        if input_text is not None:
            parameters = process.create_group("parameters")
            parameters.attrs["NX_class"] = "NXparameters"
            parameters.create_dataset(
                "sprkkr_input", data=input_text
            )
            if resolved_input_file is not None:
                parameters.create_dataset(
                    "source_file", data=str(resolved_input_file)
                )

        data = entry.create_group("data")
        data.attrs.update(NX_class="NXdata", signal="data")
        data.attrs["sprkkr_mode"] = mode
        intensity = data.create_dataset("data", data=arrays["total"])
        intensity.attrs["long_name"] = "Calculated total photoemission intensity"
        intensity.attrs["units"] = "arbitrary_units"
        if mode == "kx_ky":
            data.attrs["axes"] = ["row", "column"]
            data.attrs["row_indices"] = 0
            data.attrs["column_indices"] = 1
            data.create_dataset("row", data=np.arange(shape[0]))
            data.create_dataset("column", data=np.arange(shape[1]))
            origin = data.create_dataset("k_origin", data=k_origin)
            basis = data.create_dataset("k_basis", data=k_basis)
            origin.attrs["units"] = "1/angstrom"
            basis.attrs["units"] = "1/angstrom"
            basis.attrs["long_name"] = (
                "Momentum increments for row (NP) and column (NT)"
            )
            data.attrs["coordinate_transform"] = (
                "[kx, ky] = k_origin + row*k_basis[0] + column*k_basis[1]"
            )
            data.attrs["fixed_energy_note"] = (
                "Constant-energy momentum map: the single energy is not recorded "
                "in the native SPC columns."
            )
        elif rectilinear:
            data.attrs["axes"] = ["energy", "theta"]
            data.attrs["energy_indices"] = 0
            data.attrs["theta_indices"] = 1
            data.create_dataset("energy", data=energy[:, 0])
            data.create_dataset("theta", data=theta[0, :])
        else:
            data.attrs["axes"] = ["energy_index", "angle_index"]
            data.attrs["energy_index_indices"] = 0
            data.attrs["angle_index_indices"] = 1
            data.create_dataset("energy_index", data=np.arange(shape[0]))
            data.create_dataset("angle_index", data=np.arange(shape[1]))
            data.create_dataset("energy", data=energy)
            data.create_dataset("theta", data=theta)
            data.attrs["energy_indices"] = [0, 1]
            data.attrs["theta_indices"] = [0, 1]
        if mode == "energy":
            data["energy"].attrs["units"] = "eV"
            data["energy"].attrs["long_name"] = "Calculated initial-state energy E - E_F"
            data["theta"].attrs["units"] = "degree"
            data["theta"].attrs["long_name"] = "Photoelectron polar emission angle (SPR-KKR convention)"
        for name in ("spin_up", "spin_down", "spin_polarization", "determinant"):
            data.create_dataset(name, data=arrays[name])
        if mode == "energy":
            data.create_dataset("k_parallel", data=arrays["k_parallel"])
        data["spin_up"].attrs["units"] = "arbitrary_units"
        data["spin_down"].attrs["units"] = "arbitrary_units"
        data["spin_polarization"].attrs["units"] = "1"
        data["spin_polarization"].attrs["long_name"] = "Native spin-channel asymmetry"
        data["spin_polarization"].attrs["formula"] = "(I_1 - I_2) / (I_1 + I_2)"
        data["spin_up"].attrs["long_name"] = "SPR-KKR spin channel 1 (labelled spin up)"
        data["spin_down"].attrs["long_name"] = "SPR-KKR spin channel 2 (labelled spin down)"
        if mode == "energy":
            data["k_parallel"].attrs["units"] = "1/angstrom"
            data.attrs["k_parallel_indices"] = [0, 1]
            data["k_parallel"].attrs["long_name"] = "Native SPR-KKR k_parallel coordinate"
            if scan_k1 is not None:
                data["k_parallel"].attrs["coordinate_definition"] = (
                    "Unsigned distance in reciprocal space from KA (for 1D K-vector scans)"
                )
            elif input_parameters is not None:
                data["k_parallel"].attrs["coordinate_definition"] = (
                    "Signed k_parallel = 0.5123 sqrt(E_kin[eV]) sin(theta) [1/angstrom]"
                )
            else:
                data["k_parallel"].attrs["coordinate_definition"] = (
                    "KPARAF output; signed angular k_parallel or offset from KA in K-vector scans"
                )
        data["determinant"].attrs["long_name"] = (
            "Squared modulus of internal complex determinant (SPR-KKR diagnostic)"
        )
        if mode == "energy" and photon_energy is not None and work_function is not None:
            kinetic_energy = data.create_dataset(
                "kinetic_energy",
                data=energy + float(photon_energy) - float(work_function),
            )
            kinetic_energy.attrs["units"] = "eV"
            kinetic_energy.attrs["long_name"] = (
                "Nominal final-state electron energy relative to vacuum "
                "(negative values indicate no propagating vacuum electron)"
            )
            kinetic_energy.attrs["formula"] = "E_kin = (E-E_F) + photon_energy - work_function"
            data.attrs["kinetic_energy_indices"] = [0, 1]
        data.attrs["auxiliary_signals"] = [
            "spin_up", "spin_down", "spin_polarization", "determinant"
        ]

        source = entry.create_group("sprkkr_spc")
        source.attrs["NX_class"] = "NXcollection"
        source.create_dataset("raw_data", data=raw_data)
        source["raw_data"].attrs["columns"] = (
            ["THETA", "ENERGY", "TOTAL", "UP", "DOWN", "POLARIZATION", "K", "DETERMINANT"]
            if mode == "energy" else
            ["KX", "KY", "TOTAL", "UP", "DOWN", "POLARIZATION", "DETERMINANT"]
        )
        source.create_dataset("nt", data=int(spc.NT()))
        source.create_dataset("np", data=int(spc.NP()))
        source.create_dataset("comment", data=str(spc.COMMENT()))
    return filename
