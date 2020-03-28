import numpy as np
from collections import OrderedDict
from ase import Atoms,Atom
from ase import io
from ase.visualize import view
from ase.units import Bohr
from ase.data import chemical_symbols
import copy
import argparse

__author__='JM'
###############################################################################################
# Some usefull objects
###############################################################################################
class AtomData(object):

    def __init__(self):
        self.x=None
        self.y=None
        self.z=None
        self.typeid=None
        self.id=None
        self.symbol=None
        self.pos=np.array([0.,0.,0.])

    def __repr__(self):
        s = "<SiteData(id={:d}, pos={},{},{},type={:d})>".format(self.id, self.pos[0],self.pos[1],self.pos[2],self.typeid)
        return s

class LayerData(object):
    def __init__(self):
        self.id=1
        self.nat=1
        self.atoms=[]
        self.layvec=np.array([0.,0.,0.])
    def __repr__(self):
        s = "<LayerData(id={:d}, number of atoms in layer={:d})".format(self.id, self.nat)
        for atoms in self.atoms:
            s+= "   Atom in Layer id={:d}, pos={},{},{},type={:d} ".format(atoms.id,atoms.pos[0],atoms.pos[1],atoms.pos[2],atoms.typeid)
        s+=">\n"
        return s
###############################################################################################
# Parser for SPRKKR potential file
###############################################################################################
def _skip_lines(fd, num):
    for ii in range(num):
        line = next(fd)
    return line

def _skip_lines_to(fd, key):
    while 1:
        try:
            line = next(fd)

        except StopIteration:
            return ''

        if key in _dec(line):
            return line
class AttrDict(dict):
    """
    A dict with the attribute access to its items.
    """

    def __getattr__(self, name):
        try:
            return self[name]

        except KeyError:
            raise AttributeError(name)

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __str__(self):
        if self.keys():
            return '\n'.join(
                [self.__class__.__name__ + ':'] +
                self.format_items()
            )

        else:
            return self.__class__.__name__

    def __repr__(self):
        return self.__class__.__name__

    def __dir__(self):
        return list(self.keys())

    def copy(self):
        return type(self)(self)

    def format_items(self):
        num = max(map(len, list(self.keys()))) + 1
        return [key.rjust(num) + ': ' + repr(val)
                for key, val in sorted(self.items())]
def _dec(val):
    if isinstance(val, bytes):
        return val.decode('utf-8')

    else:
        return val

def _next_key_val(fd, dtype=bytes):
    aux = next(fd).split(maxsplit=1)
    return _dec(aux[0].lower()), dtype(aux[1].strip())

def _gen_next_key_val(fd, dtype=bytes):
    for line in fd:
        aux = line.split(maxsplit=1)
        if len(aux) == 2:
            yield _dec(aux[0].lower()), dtype(aux[1].strip())

        else:
            return

def _convert_array(x, dtype):
    return np.array([dtype(ii) for ii in _dec(x).split()])

_scf_info_dtypes = {
    int
    : {'ibzint', 'nktab', 'scf-iter'},
    float
    : {'scf-mix', 'scf-tol', 'rmsavv', 'rmsavb', 'ef', 'vmtz'},
    lambda x: _convert_array(x, int)
    : {'ne'},
}

_lattice_dtypes = {
    lambda x: [fun(x.split()[ii]) for ii, fun in enumerate([int] + 4 * [str])]
    : {'bravais'},
    float
    : {'alat'},
    lambda x: _convert_array(x, float)
    : {'a(1)', 'a(2)', 'a(3)'},
}

_sites_dtypes = {
    lambda x: _convert_array(x, float)
    : {'basscale'},
}

_reference_dtypes = {
    int
    : {'nref'},
}

_magdir_dtypes = {
    int
    : {'kmrot'},
    lambda x: _convert_array(x, float)
    : {'qmvec'},
}

def _convert(val, key, dtypes):
    for dtype, keys in dtypes.items():
        if key in keys:
            return dtype(val)

    else:
        return _dec(val)


def read_potential(filename):
    class PotFromFile(AttrDict):
        pass

    out = PotFromFile(filename=filename)

    with open(filename, 'rb') as fd:
        _skip_lines(fd, 1)
        setattr(out, *_next_key_val(fd))
        assert(out.get('header') is not None)

        _skip_lines(fd, 1)
        for _ in range(4):
            setattr(out, *_next_key_val(fd))
        assert(out.get('title') is not None)

        line = _skip_lines(fd, 2)
        assert('GLOBAL SYSTEM PARAMETER' in _dec(line))
        for _ in range(5):
            setattr(out, *_next_key_val(fd, dtype=int))

        line = _skip_lines(fd, 2)
        assert('SCF-INFO' in _dec(line))
        for key, val in _gen_next_key_val(fd):
            out[key] = _convert(val, key, _scf_info_dtypes)

        line = _skip_lines(fd, 1)
        assert('LATTICE' in _dec(line))
        for key, val in _gen_next_key_val(fd):
            out[key] = _convert(val, key, _lattice_dtypes)

        line = _skip_lines(fd, 1)
        assert('SITES' in _dec(line))
        for key, val in _gen_next_key_val(fd):
            if key == 'iq':
                out.qbas = np.loadtxt(fd, usecols=(1, 2, 3), max_rows=out.nq,
                                      ndmin=2)
                break
            out[key] = _convert(val, key, _sites_dtypes)

        line = _skip_lines(fd, 2)
        assert('OCCUPATION' in _dec(line))
        for key, val in _gen_next_key_val(fd):
            if key == 'iq':
                out.occupation = []
                for iq in range(out.nq):
                    line = next(fd).split()
                    aux1 = list(map(int, line[1:4]))
                    aux2 = [fun(line[4+ii])
                            for ii, fun in enumerate([int, float] * aux1[2])]
                    out.occupation.append(aux1 + [aux2])
                break

        line = _skip_lines(fd, 2)
        assert('REFERENCE SYSTEM' in _dec(line))
        for key, val in _gen_next_key_val(fd):
            if key == 'iref':
                out.ref = np.loadtxt(fd, usecols=(1, 2), max_rows=out.nref,
                                     ndmin=2)
                break
            out[key] = _convert(val, key, _reference_dtypes)

        line = _skip_lines(fd, 2)
        assert('HOST MADELUNG POTENTIAL' in _dec(line))
        line = _skip_lines_to(fd, 'NLMTOP-POT').split()
        out.nlmtop_pot_num = int(line[1])
        out.nlmtop_pot = np.loadtxt(
            fd, dtype=np.dtype('i4, i4, f8'),
            max_rows=out.nq * out.nlmtop_pot_num,
            ndmin=1,
        )

        line = _skip_lines(fd, 2)
        assert('CHARGE MOMENTS' in _dec(line))
        line = _skip_lines_to(fd, 'NLMTOP-CHR').split()
        out.nlmtop_chr_num = int(line[1])
        out.nlmtop_chr = np.loadtxt(
            fd, dtype=np.dtype('i4, i4, f8'),
            max_rows=out.nq * out.nlmtop_chr_num,
            ndmin=1,
        )

        line = _skip_lines(fd, 2)
        assert('MAGNETISATION DIRECTION' in _dec(line))
        # Ignores M*_T entries.
        for key, val in _gen_next_key_val(fd):
            if key == 'iq':
                aux = np.loadtxt(fd, max_rows=out.nq, ndmin=2)
                out.mag_dir = aux[:, 1:]
                break
            out[key] = _convert(val, key, _magdir_dtypes)

        _skip_lines_to(fd, 'MESH INFORMATION')
        setattr(out, *_next_key_val(fd))
        _skip_lines(fd, 1)
        out.mesh = np.loadtxt(
            fd, dtype=np.dtype([('im', 'i4'), ('r_1', 'f8'), ('dx', 'f8'),
                                ('jrmt', 'i4'), ('rmt', 'f8'), ('jrws', 'i4'),
                                ('rws', 'f8')]),
            max_rows=out.nm,
            ndmin=1,
        )

        line = _skip_lines(fd, 2)
        assert('TYPES' in _dec(line))
        line = _skip_lines(fd, 1).split()
        out.types = np.loadtxt(
            fd, dtype=np.dtype([('it', 'i4'), ('txt_t', 'U10'), ('zt', 'i4'),
                                ('ncort', 'i4'), ('nvalt', 'i4'),
                                ('nsemcorshlt', 'i4')]),
            max_rows=out.nt,
            ndmin=1,
        )

        line = _skip_lines(fd, 2)
        assert('POTENTIAL' in _dec(line))
        out.potentials_vt = []
        out.potentials_bt = []
        for ii in range(out.nt):
            line = _skip_lines(fd, 1).split()
            assert((_dec(line[0]) == 'TYPE') and
                   (int(line[1]) == out.types[ii]['it']))
            val = np.fromfile(fd, count=out.mesh[0]['jrws'], sep=' ',
                              dtype=np.float64)
            out.potentials_vt.append(val)
            val = np.fromfile(fd, count=out.mesh[0]['jrws'], sep=' ',
                              dtype=np.float64)
            out.potentials_bt.append(val)
            _skip_lines_to(fd, '========')

        line = _skip_lines(fd, 2)
        assert('CHARGE' in _dec(line))
        out.charges = []
        out.charges_bt = []
        for ii in range(out.nt):
            line = _skip_lines(fd, 1).split()
            assert((_dec(line[0]) == 'TYPE') and
                   (int(line[1]) == out.types[ii]['it']))
            val = np.fromfile(fd, count=out.mesh[0]['jrws'], sep=' ',
                              dtype=np.float64)
            out.charges.append(val)
            val = np.fromfile(fd, count=out.mesh[0]['jrws'], sep=' ',
                              dtype=np.float64)
            out.charges_bt.append(val)
            _skip_lines_to(fd, '========')

        line = _dec(_skip_lines(fd, 2)).split()
        assert('MOMENTS' in line[0])
        num = len(line[1:])
        out.moments_names = line[1:]
        out.moments = []
        for ii in range(out.nt):
            line = _skip_lines(fd, 1).split()
            assert((_dec(line[0]) == 'TYPE') and
                   (int(line[1]) == out.types[ii]['it']))
            val = np.fromfile(fd, count=num, sep=' ', dtype=np.float64)
            out.moments.append(val)
            _skip_lines_to(fd, '========')

    return out
###############################################################################################
def floatjm(inp):
        try:
                result=float(inp)
        except ValueError:
                print ("Wrong float {}, set to 0.0!".format(inp))
                result=0.0
        return result

#
###############################################################################################
ase_vis=False
# Some defaults and initialisations
filename=None
layers={}
latvec=np.zeros(shape=(2,2), dtype=float)
#READ IN COMAND LINE PARAMETERS
description='This is a sctipt to visualise in_struct.inp files. '
description+='     You can either use ase visualisation tools or export it to cif file'
parser = argparse.ArgumentParser(description=description)
parser.add_argument('-i','--input', help='in_struct.inp file name (if not specified only pot file will be ploted)',required=False)
parser.add_argument('-p','--pot', help='Filemane of scf Potential (always needed to be specified)',required=True)

parser.add_argument('-o','--out',type=str,help='Output file name for structure (default strcuture.cif)', default='strcuture.cif',required=False)
parser.add_argument('-f','--format',type=str,help='Output file fomrmat for structure (see ase allowed formats, default cif)', default='cif',required=False)
parser.add_argument('-a','--ase',help='Use ase visualisation', action='store_true',required=False)
parser.add_argument('-v','--vac',type=float,help='Size of added vacuum in AA (default=10.0)',default=10.0,required=False)
parser.add_argument('-b','--nbulk',type=int,help='Repetition of bulk unit (default=2)',default=2,required=False)


args = parser.parse_args()
Vac=args.vac
ciffile=args.out
cifpotfile=ciffile+'_pot'
outformat=args.format
ase_vis=args.ase
nbulk=args.nbulk #how many bulk repetirions will be used
filename=args.input
vis_in_structure=True
if filename == None:
   vis_in_structure=False
potfile=args.pot

Potfile=read_potential(potfile)
Pot_Atoms=Atoms()

typeid_to_symbol={}

for iq in range(Potfile.nq):
    it_iq=Potfile.occupation[iq][3][0]
    Z=Potfile.types[it_iq-1][2]
    typeid_to_symbol[iq+1]=chemical_symbols[Z]

Pot_cell=np.zeros(shape=(3,3), dtype=float)
Pot_alat=Potfile.alat*Bohr
Pot_cell[0,:]=Potfile['a(1)']
Pot_cell[1,:]=Potfile['a(2)']
Pot_cell[2,:]=Potfile['a(3)']
Pot_cell*=Pot_alat
Pot_Atoms.set_cell(Pot_cell)
Pot_pos = np.empty((0,3), float)
for iq in range(Potfile.nq):
    it_iq=Potfile.occupation[iq][3][0]
    Z=Potfile.types[it_iq-1][2]
    atm=Atom(tag=it_iq)
    atm.symbol=chemical_symbols[Potfile.types[it_iq-1][2]]
    Pot_Atoms.append(atm)
Pot_Atoms.set_positions(Potfile.qbas*Pot_alat)

if (ase_vis):
    view(Pot_Atoms)
Pot_Atoms.write(cifpotfile, format = outformat)

#From here on visualisation of in_structure.inp
if not vis_in_structure:
   exit()
# Extract data from structure file
with open(filename) as f:
    data = f.readlines()

alat=floatjm(data[1])*Bohr
nlayer=int(data[4])

for k in range(2,4):
    i=0
    for pos  in data[k].split():
        latvec[k-2,i]=floatjm(pos)
        i=i+1
line=5

for ilayer in range(nlayer):
    layer=LayerData()
    layer.nat=int(data[line])
    layer.id=ilayer+1
    for iat in range(layer.nat):
        atom=AtomData()
        line=line+1
        atom.id=int(data[line].split()[0])
        atom.typeid=int(data[line].split()[1])
        line=line+1
        x=floatjm(data[line].split()[1])
        y=floatjm(data[line].split()[2])
        z=floatjm(data[line].split()[0])
        atom.pos[0]=x
        atom.pos[1]=y
        atom.pos[2]=z
        layer.atoms.append(atom)
    line=line+1
    layers[ilayer]=layer

nvec=int(data[line].split()[0])
bulk_rep_unit=int(data[line].split()[1])
line=line+2
for ivec in range(nvec):
    line=line+1
    layers[ivec].layvec[0]=floatjm(data[line].split()[1])
    layers[ivec].layvec[1]=floatjm(data[line].split()[2])
    layers[ivec].layvec[2]=floatjm(data[line].split()[0])
    line=line+1

#
# expand list of layers in order to visualise bulk region
new_nlayer=nlayer
for i in range(nbulk):
    for k in range(bulk_rep_unit-1,nlayer):
        layers[new_nlayer]=copy.deepcopy(layers[k])
        new_nlayer+=1
# move origins of all atoms back wrt. to first layer
zmax=-99999.
zmin= 10000.
for ilay in range(1,new_nlayer):
        for atoms in layers[ilay].atoms:
                for jlay in range(0,ilay):
                        atoms.pos=atoms.pos+layers[jlay].layvec
                zmax=max(atoms.pos[2],zmax)
                zmin=min(atoms.pos[2],zmin)

# Create Atoms and Atom ASE-objects from structural data
structure=Atoms()
cell=np.zeros(shape=(3,3), dtype=float)
cell[0:2,0:2]=latvec[0:2,0:2]
cell[2,2]= zmax
cell=cell*alat
cell[2,2]= cell[2,2]+Vac #add XX AA of vacuum region
structure.set_cell(cell)
structure.set_pbc([True,True,True])
#Add atom into the Atoms
allpos = np.empty((0,3), float)
for ilay in range(new_nlayer):
    for atoms in layers[ilay].atoms:
        pos=atoms.pos.reshape(1,3)
        allpos=np.append(allpos,pos,axis=0)
        atm=Atom()
        atm.symbol=typeid_to_symbol[atoms.typeid]
        structure.append(atm)
#
# z-components needs to be transformed to the cryst. coordinates
# all others are already in cryst. cooridnates
# zmin is minum z-position and is used to shift the origin in order to avoid
# that this minimum atom will end up in the next unit cell
#
allpos[:,2]=(allpos[:,2]-zmax-zmin)*alat/cell[2,2]
structure.set_scaled_positions(allpos)
structure.set_positions(structure.get_positions(wrap=True))
# Visualise the structure or write out to the file

structure.write(ciffile, format = outformat)

if (ase_vis):
    view(structure)
