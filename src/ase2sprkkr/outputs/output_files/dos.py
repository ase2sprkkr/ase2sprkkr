from ..output_files_definitions import OutputFileValueDefinition as V, create_output_file_definition, OutputFileDefinition
from typing import Optional
import numpy as np

from ..output_files import CommonOutputFile
from ...common.grammar_types  import unsigned, Array, Table, RestOfTheFile, NumpyArray
from ...common.generated_configuration_definitions import \
    NumpyViewDefinition as NV, \
    GeneratedValueDefinition as GV
from ...visualise.plot import PlotInfo, combined_colormap, Multiplot


class DOS:
    def __init__(self, id, type, spin, l, dos):
        """ Object, that holds DOS for one atom 
        (given type and corresponding id id DOS file)
        spin, l (s,p,d,f...) """
        self.id = id
        self.type = type
        self.spin = spin
        self.l = l
        self.dos = dos

class DOSOutputFile(CommonOutputFile):

    def plot(self, layout=(2,None), figsize=(6,4), latex=True,
             filename:Optional[str]=None, show:Optional[bool]=None, dpi=600,
             **kwargs
             ):
        mp=Multiplot(layout=layout, figsize=figsize, latex=latex)
        mp.plot(self.TOTAL, **kwargs)
        mp.plot(self.UP, **kwargs)
        mp.plot(self.DOWN, **kwargs)
        mp.plot(self.POLARIZATION, **kwargs)
        mp.finish(filename, show, dpi)

    def iterate_dos(self):
        for id, t in self.TYPES():
            id =1

    def dos_for_site_type_index(self, atom):
        """ Return slice to the DOS array selecting the datas for a given site type """
        if isinstance(atom, str):
            atom=self.site_type_index(atom)
        spins = self.n_spins()
        c=0
        for i,t in enumerate(self.TYPES()):
            s=c
            c+=self.ORBITALS[t['IQAT'][0]]*spins
            if i==atom:
               return slice(s,c)

    def dos_for_site_type(self, atom, spin=None, l=None):
        """ Return density of states for a given atom,
        indexed either by integer index, or a string type.
        
        The resulting array is indexed by: [spin, l, energy], however,
        it can be restricted to given spin and/or l by arguments.
        """
        if isinstance(atom, str):
            atom=self.site_type_index(atom)
        key = self.dos_for_site_type_index(atom)
        out = self.DOS[key]
        out = out.reshape((-1, self.ORBITALS[self.TYPES[atom]['IQAT'][0]], out.shape[1]))
        out = np.moveaxis(out, 1, 0)
        if spin is not None:
           out=out[:,spin]
        if l is not None:
           out=out[l] 
        return out

    def _check_arithmetic(self, other):
        try:
            assert np.allclose(other.ENERGY(), self.ENERGY())
            assert np.allclose(other.THETA(), self.THETA())
        except AssertionError:
            raise ValueError("The outputs are not compatible to summed/subtracted.");

    def n_spins(self):
        return len(self.DOS()) // sum( (t['IQAT'][0] for t in self.TYPES()) )

class DOSDefinition(OutputFileDefinition):

    result_class = DOSOutputFile


def create_definition():

    pi = PlotInfo(
        axes = (('K', r'$k_{\parallel} $(${\rm \AA}^{-1}$)'),
                ('ENERGY', r'$E-E_{\rm F}$ (eV)')),
        #axes = ('K', 'ENERGY' ),
        show_zero_line = True,
        colormap = combined_colormap()
    )

    def i(j):
        return slice(None),j



    definition = create_output_file_definition('DOS', [
      V('DOS-FMT', str),
      V('RAW_DATA', NumpyArray(written_shape=(-1,8),
                               delimiter=10,indented=(80,10),
                               item_format='%8.4E'
                    ), name_in_grammar=False),

      NV('ENERGY', 'RAW_DATA', i(0), info='Energies'),
      NV('Y', 'RAW_DATA', i(1), info='Y'),
      NV('DOS', 'RAW_DATA', i(slice(2,None)), ('NE', -1),
                transpose=True,
                transform_key= lambda k,dos: dos.dos_for_site_type_index(k)\
                                             if isinstance(k, str) else k,
                info='Desntity of states', description=
                "Density of states, the leading dimension iterates: "
                "1..n_atoms, 1..n_spins, 1..n_orbitals(atoms) "
    ),
    ], cls=DOSDefinition)
    return definition

definition = create_definition()
