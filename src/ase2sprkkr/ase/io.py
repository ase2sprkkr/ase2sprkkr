from .. import Potential


def read_sprkkr(fd):
    return Potential.from_file(fd).atoms


def write_sprkkr(fd, atoms):
    pot = Potential.from_atoms(atoms)
    return pot.save_to_file(fd)
