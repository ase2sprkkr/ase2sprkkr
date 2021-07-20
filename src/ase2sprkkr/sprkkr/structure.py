from ase import Atom
from .sprkkr_atoms import SprKkrAtoms
from ..potential.potentials import Potential
from ase.units import Bohr
import copy
from ase.data import chemical_symbols

"""
This file contains function structure_file_to_atoms, which is used for
visualisation of a surface structure, using .pot and in_structure.inp
files.

TODO: This implementaion can handle only one-purpose reading of the structure
file.
"""

### Helper objects and functions

class AtomData(object):
    """ A helper object for reading in_structure.inp file """
    def __init__(self):
        self.x=None
        self.y=None
        self.z=None
        self.typeid=None
        self.id=None
        self.symbol=None
        self.pos=np.array([0.,0.,0.])

    @staticmethod
    def from_text(lines):
        """ Create the Atom data object from the two text lines
        :from a in_structure.inp file """
        atom = AtomData()
        atom.id, atom.type_id = map(int, data[0].split()[:2])
        atom.pos[:] = map(floatjm, data[1].split()[:3])

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

  atoms: ase2sprkkr.sprkkr.sprkkr_atoms.SprKkrAtoms
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
          atom=AtomData.from_text(data[line:line+2])
          line+=2
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
  for i in range(n_bulk):
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
  structure=SprKkrAtoms()
  cell=np.zeros(shape=(3,3), dtype=float)
  cell[0:2,0:2]=latvec[0:2,0:2]
  cell[2,2]= zmax
  cell=cell*alat
  cell[2,2]= cell[2,2]+vacuum_height #add XX AA of vacuum region
  structure.set_cell(cell)
  structure.set_pbc([True,True,True])
  #Add atom into the Atoms
  num_atoms = 0
  for layer in layers:
      num_atoms += len(layer.atoms)
  allpos = np.empty((num_atoms,3), float)
  i = 0

  for layer in layers:
      for atom in layer.atoms:
          allpos[i,:] = atom.pos
          atm=Atom()
          atm.symbol= chemical_symbols[ potential.TYPES.DATA()[atom.typeid - 1]['ZT'] ]
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

__all__ = [ structure_file_to_atoms ]
