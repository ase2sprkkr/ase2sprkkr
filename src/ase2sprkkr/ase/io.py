from .. import Potential


def read_sprkkr(fd, index=None, **kwargs):
    out = Potential.from_file(fd).atoms
    if index is not None:
        return [out][index]
    return out


read_SPRKKR = read_sprkkr


def write_sprkkr(fd, atoms, **kwargs):
    pot = Potential.from_atoms(atoms[0])
    return pot.save_to_file(fd)


write_SPRKKR = write_sprkkr
