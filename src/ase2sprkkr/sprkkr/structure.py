"""
Helper classes for a2s_visualise_in_struct utility.

This file contains function structure_file_to_atoms, which is used for
visualisation of a surface structure, using .pot and in_structure.inp
files.

TODO: This implementaion can handle only one-purpose reading of the structure
file.
"""


from ase import Atom
from .sprkkr_atoms import SPRKKRAtoms
from ..potentials.potentials import Potential
from ase.units import Bohr
import copy
from ase.data import chemical_symbols
import numpy as np

### Helper objects and functions

class AtomData(object):
    """ A helper object for reading in_structure.inp file """
    def __init__(self):
        self.x=None
        self.y=None
        self.z=None
        self.type_id=None
        self.id=None
        self.symbol=None
        self.pos=np.array([0.,0.,0.])

    @staticmethod
    def from_text(lines):
        """ Create the Atom data object from the two text lines
        :from a in_structure.inp file """
        atom = AtomData()
        atom.id, atom.type_id = map(int, lines[0].split()[:2])
        atom.pos[:] = list(map(floatjm, lines[1].split()[:3]))
        return atom

    def __repr__(self):
        s = "<SiteData(id={:d}, pos={},{},{},type={:d})>".format(self.id, self.pos[0],self.pos[1],self.pos[2],self.type_id)
        return s

    def get_symbol(self, potential):
        IT = potential.OCCUPATION.DATA()[self.type_id - 1][3][0][0]
        z = potential.TYPES.DATA()[IT - 1]['ZT']
        return chemical_symbols[z]
class LayerData(object):
    def __init__(self):
        self.id=1
        self.nat=1
        self.atoms=[]
        self.layvec=np.array([0.,0.,0.])

    def __repr__(self):
        s = "<LayerData(id={:d}, number of atoms in layer={:d})".format(self.id, self.nat)
        for atom in self.atoms:
            s+= "   Atom in Layer id={:d}, pos={},{},{},type={:d} ".format(atom.id,atom.pos[0],atom.pos[1],atom.pos[2],atom.type_id)
        s+=">\n"
        return s

def floatjm(inp):
        try:
                result=float(inp)
        except ValueError:
                print ("Wrong float {}, set to 0.0!".format(inp))
                result=0.0
        return result

###############################################################################################


def structure_file_to_atoms(filename, potential:Potential, n_bulk:int=2, vacuum_height:float=10.0):
  """
  Read in_structure.inp file (that contains informations about the structure of a surface) 
  and create the ASE atoms object according to the readed data.

  Parameters
  ----------
  filename: str
    File to read

  atoms: ase2sprkkr.sprkkr.sprkkr_atoms.SPRKKRAtoms
    Atoms object, from where the types of atoms will be given

  n_bulk: int
    Number of repetition of bulk atoms

  vacuum_height: float
    Height of the vacuum above the surface.
  """
  with open(filename) as f:
      data = f.readlines()

  alat=floatjm(data[1])*Bohr
  nlayer=int(data[4])

  latvec=np.zeros(shape=(2,2), dtype=float)
  for k in range(2,4):
      i=0
      for pos in data[k].split():
          latvec[k-2,i]=floatjm(pos)
          i=i+1
  line=5

  layers={}
  for ilayer in range(nlayer):
      layer=LayerData()
      layer.nat=int(data[line])
      layer.id=ilayer+1
      for iat in range(layer.nat):
          line+=3
          atom=AtomData.from_text(data[line-2:line])
          layer.atoms.append(atom)
      layers[ilayer]=layer

  nvec=int(data[line].split()[0])
  bulk_rep_unit=int(data[line].split()[1])
  line=line+2
  for ivec in range(nvec):
      line+=1
      layers[ivec].layvec[0]=floatjm(data[line].split()[1])
      layers[ivec].layvec[1]=floatjm(data[line].split()[2])
      layers[ivec].layvec[2]=floatjm(data[line].split()[0])
      line+=1

  #
  # expand list of layers in order to visualise bulk region
  new_nlayer=nlayer
  for i in range(n_bulk):
      for k in range(bulk_rep_unit-1,nlayer):
          layers[new_nlayer]=copy.deepcopy(layers[k])
          new_nlayer+=1
  # move origins of all atoms back wrt. to first layer
  zmax=-99999.
  zmin= 10000.
  for ilay in range(1,new_nlayer):
          for atom in layers[ilay].atoms:
                  for jlay in range(0,ilay):
                          atom.pos=atom.pos+layers[jlay].layvec
                  zmax=max(atom.pos[2],zmax)
                  zmin=min(atom.pos[2],zmin)

  # Create Atoms and Atom ASE-objects from structural data
  structure=SPRKKRAtoms()
  cell=np.zeros(shape=(3,3), dtype=float)
  cell[0:2,0:2]=latvec[0:2,0:2]
  cell[2,2]= zmax
  cell=cell*alat
  cell[2,2]= cell[2,2]+vacuum_height #add XX AA of vacuum region
  structure.set_cell(cell)
  structure.set_pbc([True,True,True])
  #Add atom into the Atoms
  num_atoms = 0
  for layer in layers.values():
      num_atoms += len(layer.atoms)
  allpos = np.empty((num_atoms,3), float)
  i = 0

  for layer in layers.values():
      for atom in layer.atoms:
          allpos[i,:] = atom.pos
          i+=1
          atm=Atom()
          atm.symbol=atom.get_symbol(potential)
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

  return structure

__all__ = [ 'structure_file_to_atoms' ]
