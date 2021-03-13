from ase import Atoms

class SprKkrAtoms(Atoms):
   """ ASE Atoms object extended by the data necessary for SPR-KKR calculations """

   @staticmethod
   def convert_from_ase_atoms(obj):
       """ Convert ASE Atoms object to the one usable by 
           SprKkr. 
       """
       if obj:
          obj.__class__ = SprKkrAtoms
       return obj
