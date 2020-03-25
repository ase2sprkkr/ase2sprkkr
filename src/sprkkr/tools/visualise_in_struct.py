import numpy as np
from collections import OrderedDict
from ase import Atoms,Atom
from ase import io
from ase.visualize import view
from ase.units import Bohr
import copy

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

filename='in_structur.inp2'

with open(filename) as f:
    data = f.readlines()
# Extract data
alat=float(data[1])*Bohr
nlayer=int(data[4])
layers={}

latvec=np.zeros(shape=(2,2), dtype=float)


for k in range(2,4):
    i=0
    for pos  in data[k].split():
        latvec[k-2,i]=float(pos)
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
        x=float(data[line].split()[1])
        y=float(data[line].split()[2])
        z=float(data[line].split()[0])
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
    layers[ivec].layvec[0]=float(data[line].split()[1])
    layers[ivec].layvec[1]=float(data[line].split()[2])
    layers[ivec].layvec[2]=float(data[line].split()[0])
    line=line+1

#
typeid_to_symbol={}
for key,lay in layers.items():
    for atm in lay.atoms:
        typeid_to_symbol[atm.typeid]=''

typeid_to_symbol[3]='Te'
typeid_to_symbol[1]='Hf'
typeid_to_symbol[2]='Te'
typeid_to_symbol[4]='H'

nbulk=5 #how many bulk repetation will be used
# expand list of layers in order to visualise bulk region
new_nlayer=nlayer
for i in range(nbulk):
    for k in range(bulk_rep_unit-1,nlayer):
        layers[new_nlayer]=copy.deepcopy(layers[k])
        new_nlayer+=1
# move origins of all atoms back wrt. to first layer
zmax=-99999.

for ilay in range(1,new_nlayer):
        for atoms in layers[ilay].atoms:
                for jlay in range(0,ilay):
                        atoms.pos=atoms.pos+layers[jlay].layvec
                zmax=max(atoms.pos[2],zmax)

# Create Atoms and Atom objects from ASE
print(zmax)
structure=Atoms()
cell=np.zeros(shape=(3,3), dtype=float)
cell[0:2,0:2]=latvec[0:2,0:2]
cell[2,2]= zmax
cell=cell*alat
cell[2,2]= cell[2,2]+10.0 #add 10 AA of vacuum region
print(cell)
structure.set_cell(cell)
structure.set_pbc([True,True,True])
#Add atom into the Atoms
allpos = np.empty((0,3), float)
for ilay in range(new_nlayer):
    for atoms in layers[ilay].atoms:
        pos=atoms.pos.reshape(1,3)
        allpos=np.append(allpos,pos,axis=0)
#        allpos=np.append(allpos,pos,axis=0)
        atm=Atom()
        atm.symbol=typeid_to_symbol[atoms.typeid]
        structure.append(atm)
#z-components needs to be transformed to the cryst. coordinates
#all others are already in cryst. cooridnates
allpos[:,2]=allpos[:,2]*alat/cell[2,2]
print(allpos[:,2])
structure.set_scaled_positions(allpos)

structure.write('str.cif', format = 'cif')

view(structure)
