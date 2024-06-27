from ..output_files_definitions import OutputFileValueDefinition as V, create_output_file_definition, OutputFileDefinition
from typing import Optional
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import copy

from ..output_files import CommonOutputFile, Arithmetic
from ...common.grammar_types import NumpyArray
from ...common.decorators import cached_property
from ...common.generated_configuration_definitions import \
    NumpyViewDefinition as NV
from ...visualise.plot import Multiplot

from ase.units import Rydberg
from packaging.version import Version


class DOS(Arithmetic):

    def __init__(self, energy, dos, type=None, id=None, spin=None, l=None):  # NOQA
        """ Object, that holds DOS for one atom
        (given type and corresponding id id DOS file)
        spin, l (s,p,d,f...) """
        self.id = id
        self.type = type
        self.energy = energy

        self.spin = spin
        self.l = l    # NOQA
        self.dos = dos

    def copy(self, copy_values=True):
        out = copy.copy(self)
        if copy_values:
            out.dos = out.dos.copy()
        return out

    @property
    def shape(self):
        return self.dos.shape

    def __getitem__(self, key):
        return self.dos[key]

    def plot(self, axis=None, legend_ncols=1, legend_height=0.3, legend_fontsize=10, **kwargs):
        axis.set_xlabel(r'$E-E_{\rm F}$ (eV)')
        axis.set_ylabel(r'DOS (states/eV)')

        title = []
        if self.type:
           title=[ self.type ]
        if self.spin is not None:
           title.append('spin {}'.format({0:'up', 1:'down'}.get(self.spin, self.spin)))
        if self.l is not None:
           orbital = {0:'s',1:'p', 2:'d', 3: 'f'}.get(self.l, None)
           title.append(orbital + 'orbitals' if orbital else self.l)
        if title:
           axis.title.set_text(', '.join(title))

        params = {
            's': { 'color': 'blue' },
            'p': { 'color': 'green' },
            'd': { 'color': 'orange' },
            'f': { 'color': 'cyan' },
            'total' : {'color': 'black' }
        }

        def plot_l(data, spin, l):  # NOQA
            i = l if l in params else 'total'
            args = params[i]
            line, = axis.plot(self.energy, data, label=l, **args)
            return line

        def plot_spin(data, spin, legend):
            if spin == -1:
               data = data * spin
            if self.l is not None:
               handles = [ plot_l(data, spin, self.l) ]
            else:
               handles = [
                   plot_l(d, spin, l) for
                   d,l in zip(data, ('s','p','d','f'))
               ]
               if len(handles):
                 handles.append(
                   plot_l(np.sum(data, axis=0), spin, 'total')
                 )
            if legend and not self.l:
                ncols = 'ncols' if Version(matplotlib.__version__) >= Version("3.6") else 'ncol'
                axis.legend(handles=handles,loc='best',fontsize=legend_fontsize, **{ncols : legend_ncols }, handleheight=legend_height)

        if self.spin is not None:
           plot_spin(self.dos, self.spin or 1, True)
        else:
           legend=True
           for d,s in zip(self.dos, (1,-1)):
               plot_spin(d,s, legend=legend)
               legend=False

    def _do_arithmetic(self, func, other):
        if isinstance(other, DOS):
            other = other.dos
        getattr(self.dos, func)(other)
        self.id = None
        self.type = None

    def _check_arithmetic(self, other):
        if self.energy is other.energy:
            return
        assert np.allclose(self.energy, other.energy)

    def __repr__(self):
        if self.type:
           return f'<DOS of {self.type}>'
        else:
           return '<DOS>'


class DOSOutputFile(CommonOutputFile):

    def __init__(self, definition, container=None):
        super().__init__(definition, container)
        self.ENERGY.add_hook(self._clear_computed)
        self.EFERMI.add_hook(self._clear_computed)

    def _clear_computed(self, _):
        if 'energy' in self.__dict__:
            del self.energy

    @cached_property
    def energy(self):
        return (self.ENERGY() - self.EFERMI()) * Rydberg

    def plot(self, spin=None, l=None, layout=2, figsize=(6,4), latex=True,  # NOQA
             filename:Optional[str]=None, show:Optional[bool]=None, dpi=600,
             **kwargs
             ):
        if isinstance(layout, int):
          layout = ( (self.n_types() ) // layout +1, layout)
        mp=Multiplot(layout=layout, figsize=figsize, latex=latex)
        plt.subplots_adjust(left=0.12,right=0.95,bottom=0.1,top=0.9, hspace=0.6, wspace=0.4)
        to_plot = self.iterate_dos(spin, l, total=True)
        for dos, axis in zip(to_plot, mp):
            dos.plot(axis)
        mp.finish(filename, show, dpi)

    def __getitem__(self, name):
        """
        In addition to the general Container __getitem__, the following
        is possible>
        ``dos[0].plot()`` - plot the first atom
        ``dos['Te'].plot()`` - plot Te atom
        """
        try:
            return super().__getitem__(name)
        except KeyError as ke:
            if not np.issubdtype(name.__class__, np.integer):
                for i, type in enumerate(self.TYPES):
                  if type[0] == name:
                      name = i
                      break
                else:
                  raise KeyError(f"There is no atomic type {name} nor such value"
                                 " in the DOS file") from ke
            for i,slc in enumerate(self.iterate_data_slices()):
                if i == name:
                     return self._create_dos(slc, i)
            raise KeyError(f"There is no {name}th atomic type"
                           " in the DOS file") from ke

    def total_dos(self, spin=None, ll=None):
        return [i for i in self.iterate_dos(spin, ll, total=True) ][-1]

    def __iter__(self):
        return self.iterate_dos(total=False)

    @staticmethod
    def _resolve_spin(spin):
        if isinstance(spin, str):
            return {'up': 0, 'down':1}[spin.lower()]
        return spin

    def iterate_data_slices(self):  # NOQA
        spins = self.n_spins()
        end=0

        for t in self.TYPES():
            start=end
            orbitals=self.n_orbitals_for(t)
            end+=orbitals * spins
            yield slice(start,end)

    def iterate_dos(self, spin=None, l=None, total=True):  # NOQA
        spin = self._resolve_spin(spin)
        total = bool(total)
        for i,slic in enumerate(self.iterate_data_slices()):
            out = self._create_dos( slic, i, spin, l)
            yield out
            if not total:
                continue
            if total is True:
                total = out * self.TYPES[i][2]
            else:
                total.dos += out.dos * self.TYPES[i][2]
        if total:
            total.type = 'total'
            yield total

    def index_of_dos_for_site_type(self, atom):
        """ Return slice to the DOS array selecting the datas for a given site type """
        if isinstance(atom, str):
            atom=self.site_type_index(atom)
        for i,slic in enumerate(self.iterate_data_slices()):
            if i==atom:
               return slic

    def dos_for_site_type(self, atom, spin=None, l=None):  # NOQA
        """ Return density of states for a given atom,
        indexed either by integer index, or a string type.

        The resulting array is indexed by: [l, spin, energy], however,
        it can be restricted to given spin and/or l by arguments.
        """
        spin = self._resolve_spin(spin)
        if isinstance(atom, str):
            atom=self.site_type_index(atom)
        key = self.index_of_dos_for_site_type(atom)
        return self._create_dos(key, atom, spin, l)

    def _create_dos(self, key, id, spin=None, l=None):  # NOQA
        type = self.TYPES[id]['TXT_T']
        out = self.DOS[key]
        out = out.reshape((-1, self.n_orbitals_for(self.TYPES[id]), out.shape[1]))
        # out = np.moveaxis(out, 1, 0)
        if l is not None:
           out=out[:,l]
        if spin is not None:
           out=out[spin]
        return DOS( self.energy, out / Rydberg, type, id, spin, l)

    def n_orbitals_for(self, type):
        """
        Return the number of orbitals for the given type record
        """
        return max(
            self.ORBITALS[iq - 1] for iq in type['IQAT']
        )

    def n_spins(self):
        """
        Return the number of spins for each orbital
        """
        ln = len(self.DOS())
        orbitals = sum( self.n_orbitals_for(t) for t in self.TYPES() )
        return ln // orbitals


class DOSDefinition(OutputFileDefinition):

    result_class = DOSOutputFile


def create_definition():

    def i(j):
        return slice(None),j

    definition = create_output_file_definition('DOS', [
      V('DOS-FMT', str, written_name='DOS-FMT:'),
      V('RAW_DATA', NumpyArray(written_shape=(-1,8),
                               delimiter=10,indented=(80,10),
                               item_format='%8.4E'
                    ), name_in_grammar=False),

      NV('ENERGY', 'RAW_DATA', i(0), info='Energies'),
      NV('Y', 'RAW_DATA', i(1), info='Y'),
      NV('DOS', 'RAW_DATA', i(slice(2,None)), ('NE', -1),
                transpose=True,
                transform_key=lambda k,dos: dos.index_of_dos_for_site_type(k) if
                                            isinstance(k, str) else k,
                info='Desntity of states', description=
                "Density of states, the leading dimension iterates: "
                "1..n_atoms, 1..n_spins, 1..n_orbitals(atoms) "
    ),
    ], cls=DOSDefinition, info='Result of a DOS (density of states) calculation.')
    return definition


definition = create_definition()
