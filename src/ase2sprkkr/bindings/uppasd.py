from enum import Enum
from pathlib import Path
from typing import Optional
import warnings
import numpy as np
from ase import Atoms
from ..sprkkr.sprkkr_atoms import SPRKKRAtoms
from ..output_files.definitions.jxc import JXCOutputFile, Coordinates

POSFILE_FORMAT = "{:5.0f} {:5.0f} {: 4.10f} {: 4.10f} {: 4.10f}\n"
MOMFILE_FORMAT = "{:5.0f} {:5.0f} {: 4.10f} {: 4.10f} {: 4.10f} {: 4.10f}\n"
JFILE_FORMAT = "{:5.0f}  {:5.0f}  {: 4.10f}  {: 4.10f}  {: 4.10f}  {: 4.8f}  {:4.10f}\n"
DMFILE_FORMAT = "{:5.0f}  {:5.0f}  {: 4.10f}  {: 4.10f}  {: 4.10f}  {: 4.8f}  {: 4.8f}  {: 4.8f}  {:4.10f}\n"
COORDINATE_FIELDS = {
    Coordinates.cartesian: ("DRX", "DRY", "DRZ"),
    Coordinates.lattice: ("N1", "N2", "N3"),
}


def write_inpsd_file(
    atoms: Atoms, file_name=None, directory=None, id=None, sym=0, ncell=None, **kwargs
):
    file_name = file_name or "inpsd.dat"
    if directory:
        file_path = Path(directory) / file_name
    else:
        file_path = Path(file_name)

    bc = ["P" if i else "0" for i in atoms.pbc]
    bc = "            ".join(bc)
    cell = ""
    for row in atoms.cell.array:
        cell += "    {:>14.10f} {:>14.10f} {:>14.10f}\n".format(*row)
    if isinstance(ncell, int):
        ncell = (ncell, ncell, ncell)
    if ncell is not None:
        ncell = f"ncell {ncell[0]:<13} {ncell[1]:<13} {ncell[2]:<13}   System size"
    else:
        ncell = ""

    def formated(v):
        if isinstance(v, (bool, np.bool)):
            return "Y" if v else "N"
        return str(v)

    with file_path.open("w") as infile:
        infile.write(
            f"""
simid {id or str(atoms.symbols)}
{ncell}
BC    {bc}  Boundary conditions (0=vacuum, P=periodic)
cell  {cell[6:]}
Sym   {sym}    Symmetry of exchange bonding vectors (0 for no, 1 for cubic, 2 for 2d cubic, 3 for hexagonal)

"""
        )
        for k, v in kwargs.items():
            infile.write(f"{k:<14} {formated(v)}\n")


def write_jfile(
    jxc_ouput_file: JXCOutputFile,
    file_name=None,
    directory=None,
    selector=None,
    iq=None,
    it=None,
    exclude_it=None,
    exclude_vc=True,
    exchange_radius=None,
    coordinates: Coordinates = Coordinates.lattice,
):
    coordinate_fields = COORDINATE_FIELDS[coordinates]
    rows = jxc_ouput_file.filtered_data(
        selector=selector,
        iq=iq,
        it=it,
        exclude_it=exclude_it,
        exclude_vc=exclude_vc,
        exchange_radius=exchange_radius,
    )

    file_name = file_name or "jfile.dat"
    if directory:
        file_path = Path(directory) / file_name
    else:
        file_path = Path(file_name)

    with file_path.open("w") as jfile:
        for row in rows:
            jfile.write(
                JFILE_FORMAT.format(
                    row["IQ"],
                    row["JQ"],
                    row[coordinate_fields[0]],
                    row[coordinate_fields[1]],
                    row[coordinate_fields[2]],
                    row["JXX"],
                    row["DR"],
                )
            )
    return True


def write_dmfile(
    jxc_ouput_file: JXCOutputFile,
    file_name="dmfile.dat",
    directory=None,
    selector=None,
    iq=None,
    it=None,
    exclude_it=None,
    exclude_vc=True,
    exchange_radius=None,
    coordinates: Coordinates = Coordinates.lattice,
):
    coordinate_fields = COORDINATE_FIELDS[coordinates]
    rows = jxc_ouput_file.filtered_data(
        selector=selector,
        iq=iq,
        it=it,
        exclude_it=exclude_it,
        exclude_vc=exclude_vc,
        exchange_radius=exchange_radius,
    )

    if directory:
        file_path = Path(directory) / file_name
    else:
        file_path = Path(file_name)

    with file_path.open("w") as dmfile:
        for row in rows:
            dmfile.write(
                DMFILE_FORMAT.format(
                    row["IQ"],
                    row["JQ"],
                    row[coordinate_fields[0]],
                    row[coordinate_fields[1]],
                    row[coordinate_fields[2]],
                    row["DX"],
                    row["DY"],
                    row["DZ"],
                    row["DR"],
                )
            )
    return True


def write_pos_file(
    atoms: Atoms,
    file_name="posfile.dat",
    directory=None,
    selector=None,
    iq=None,
    it=None,
    exclude_it=None,
    exclude_vc=True,
    jxc_ouput_file: Optional[JXCOutputFile] = None,
):
    atoms = SPRKKRAtoms.promote_ase_atoms(atoms)
    if not jxc_ouput_file:
        jxc_ouput_file = JXCOutputFile.from_atoms(atoms)
    positions = atoms.get_scaled_positions(wrap=False)
    selected_iqs = jxc_ouput_file.iq_selector(
        selector=selector, iq=iq, it=it, exclude_it=exclude_it, exclude_vc=exclude_vc
    )
    if selected_iqs is ...:
        selected_iqs = range(1, len(atoms.sites) + 1)
    elif not selected_iqs:
        return False
    site_types = {}

    if directory:
        file_path = Path(directory) / file_name
    else:
        file_path = Path(file_name)

    with file_path.open("w") as posfile:
        for iq in selected_iqs:
            index = iq - 1
            site = atoms.sites[index]
            position = positions[index]
            site_type = site.site_type
            try:
                site_kind = site_types[site_type]
            except KeyError:
                site_kind = site_types[site_type] = len(site_types) + 1
            posfile.write(
                POSFILE_FORMAT.format(
                    iq,
                    site_kind,
                    position[0],
                    position[1],
                    position[2],
                )
            )
    return True


def write_mom_file(
    atoms: Atoms,
    file_name="momfile.dat",
    directory=None,
    selector=None,
    iq=None,
    it=None,
    exclude_it=None,
    exclude_vc=True,
    jxc_ouput_file: Optional[JXCOutputFile] = None,
):
    atoms = SPRKKRAtoms.promote_ase_atoms(atoms)
    if not jxc_ouput_file:
        jxc_ouput_file = JXCOutputFile.from_atoms(atoms)

    site_types = []
    out = []
    warned = False

    selector = jxc_ouput_file.create_selector(
        selector=selector, iq=iq, it=it, exclude_it=exclude_it, exclude_vc=exclude_vc
    )
    iq_selector = selector.iq
    it_selector = selector.it

    # moments are written in the order of atomic types
    out = {}
    it = 0

    for id, site in enumerate(atoms.sites):
        iq = id + 1
        if iq_selector is not ... and iq not in iq_selector:
            continue
        for at in site.occupation:
            if at.moments is None:
                if not warned:
                    warnings.warn(
                        f"Some atomic types do not have moments computed. Skipping."
                    )
                    warned = True
                continue
            iqs = out.get(at, None)
            if iqs is False:
                continue
            if iqs is None:
                it += 1
                if it_selector is not ... and it not in it_selector:
                    out[at] = False
                    continue
                out[at] = iqs = [iq]
            else:
                iqs.append(iq)

    if not out:
        return False

    if directory:
        file_path = Path(directory) / file_name
    else:
        file_path = Path(file_name)

    with file_path.open("w") as momfile:
        for atomic_type, iqs in out.items():
            moment = atomic_type.moments.spin_moment
            for iq in iqs:
                momfile.write(MOMFILE_FORMAT.format(iq, 1, moment, 0, 0, 1))

    return True
