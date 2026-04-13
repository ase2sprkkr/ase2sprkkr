from ..output_files_definitions import OutputFileValueDefinition as V, create_output_file_definition, OutputFileDefinition,Separator, BlankSeparator, OutputFileSectionDefinition
from typing import Optional
import numpy as np
import matplotlib.colors as mcolors
from functools import lru_cache
import warnings
from ase.units import Rydberg
from ase2sprkkr.physics.broadening import create_lorentz_broadener, create_gaussian_broadener, compute_wlortab
from scipy.interpolate import interp1d

from ..output_files import Arithmetic, CommonOutputFile
from ...common.grammar_types import Array, NumpyArray, Keyword, Char, Table, \
                                    Sequence, Complex, Real
from ...common.generated_configuration_definitions import GeneratedValueDefinition
from ...common.configuration_definitions import KeywordSeparator, SeparatorDefinition, BaseDefinition
from ...gui.plot import Multiplot, plotting_function, set_up_common_plot, single_plot
from ...common.decorators import cached_property, add_to_signature

class RATOutputFile(CommonOutputFile):

    extension = 'rat'
    plot_parameters = {}

    """
    Output file for Bloch spectral functions
    """

    def energies(self, group=0):
        """ Energies in Ev relative to E-Fermi."""
        efermi = self.EFERMI()
        return np.array([(i.ENERGY()[0].real - efermi) * Rydberg for i in self.GROUPS[group].DATA.values()])

    def order(self):
        """ Order of the energie dataset."""

    def validate_data(self):
        energies = self.energies()
        for i in range(1, len(self.GROUPS)):
            if not(self.energies(i) != energies):
                raise ValueError(f"Set of energies for GROUP={i} differs from the ones for GROUP=0")
        cst = self.core_state_types
        if len(cst[0]) > 2:
            raise ValueError("Too much KAP core-state-types")
        if len(cst[0])==2:
            itr = iter(cst[1])
            changes=0
            prev = next(i)
            for i in itr:
                if prev!=i:
                    changes+=1
            if changes >= len(cst[0]):
                raise ValueError("KAP types can't be interleaved")

    @cached_property
    def core_states_index(self):
        """ Return indexes of core_states in CORE_STATES table, some of them can have two lines """
        STATES = self.GROUPS[0].CORE_STATES
        return [ i for i in range(len(STATES)) if i == 0 or STATES[i]['ICST'] != STATES[i-1]['ICST'] ]

    @cached_property
    def core_state_types(self):
        """ classify the core states according to the KAP property """
        return np.unique(self.GROUPS[0].CORE_STATES['KAP'], return_inverse=True)

    def core_states_energies(self, no_core_splitting=False):
        states_index = self.core_states_index
        n_groups = len(self.GROUPS)
        efermi = self.EFERMI()
        out = np.empty((n_groups, len(states_index)))

        def fill_group(i):
            out[i] = self.GROUPS[i].CORE_STATES['E(Ry)'][states_index] * Rydberg

        if no_core_splitting:
            fill_group(0)
            out[1:] = out[0]
        else:
            for i in range(n_groups):
                fill_group(i)

        return out

    @cached_property
    def mendeleev(self):
        import mendeleev
        return mendeleev.element(self.TYPES[self.GROUPS[0].IT()-1]['TXT_T'])

    def lorentz_width(self, source='campbell-papp'):
        """
        Return tuple of Lorentz widths (for the two possible KAP types).
        The `source` argument is forwarded to the core_hole_width lookup.
        """
        from ase2sprkkr.physics.core_hole_width import core_hole_width
        atomic_number = self.mendeleev.atomic_number
        nc = self.GROUPS[0].NCXRAY()
        lc = self.GROUPS[0].LCXRAY()
        return (
            core_hole_width(atomic_number, nc, lc, 1, source=source),
            core_hole_width(atomic_number, nc, lc, 2, source=source),
        )

    def broadener(self, energies, gauss_width, lorentz_width, n_valence=None, wlortab=None):
        if lorentz_width is None:
            lorentz_width = self.lorentz_width()
        if n_valence is None:
            n_valence = self.mendeleev.nvalence()
        # Accept Python floats and NumPy scalar floats here
        if np.isscalar(lorentz_width):
            lorentz_width = [ float(lorentz_width) ]

        # Allow user to pass a precomputed wlortab array; otherwise compute from n_valence
        if wlortab is None:
            wlortab = compute_wlortab(energies, n_valence)

        lorentz = [
            None if w < 0.001 and n_valence > 0 else
            create_lorentz_broadener(energies, w, wlortab) for w in lorentz_width
        ]
        if len(lorentz) == 1:
            lorentz*=2

        if gauss_width > 0.001:
            gauss = create_gaussian_broadener(energies, gauss_width, method='dense')
        else:
            gauss = None

        def broaden(data, k):
            if lorentz[k]:
                data=lorentz[k](data)
            if gauss:
                data=gauss(data)
            return data

        return broaden

    @property
    def data(self):
        if not hasattr(self, '_data'):
            raise ValueError('First call the `generate_data` method to compute the data')
        return self._data

    def plot(self, layout=(None,2), figsize=None, latex=None,
             filename:Optional[str]=None, show:Optional[bool]=None, dpi=800,
             separate_plots=False,
             interpolate_to_fermi=True, interpolation_threshold=1e-4,
             zero_below=True, zero_below_num=50, zero_below_energy=25,
             no_core_splitting=False, merge_all=None, updown_layout=True,
             gauss_width=0.1, lorentz_width=None, n_valence=None,
             core_hole_width='campbell-papp',
             **kwargs):

       data = self.generate_data(interpolate_to_fermi, interpolation_threshold,
             zero_below, zero_below_num, zero_below_energy,
             no_core_splitting, merge_all,
             gauss_width, lorentz_width, n_valence, core_hole_width)

       num = data['POLARIZATION'].shape[0];   #1 or 2
       if 'SPIN' in data:
            num+=1
       num=num*2-1
       if layout[0] is None:
          layout = ( num // layout[1] +1, layout[1] )
       elif layout[1] is None:
          layout = ( layout[0],  num // layout[0] )

       if figsize is None:
          figsize = (layout[0] * 4, layout[1] * 4 )

       with Multiplot(layout=layout, figsize=figsize, latex=latex,
                       filename=filename, show=show, dpi=dpi, updown_layout=updown_layout, separate_plots=separate_plots,
                       adjust={'left':0.12, 'right':0.95, 'bottom':0.17, 'top':0.90, 'hspace':0.75, 'wspace':0.5},
                       **kwargs) as mp:
            self.POLARIZATION.plot(_inside_plot=mp,  **kwargs)
            self.DIFFERENCE.plot(_inside_plot=mp, **kwargs)
            if 'SPIN' in data:
                self.SPIN.plot(_inside_plot=mp, **kwargs)
                self.ORBIT.plot(_inside_plot=mp, **kwargs)

    def generate_data(self, interpolate_to_fermi=True, interpolation_threshold=1e-4,
             zero_below=True, zero_below_num=50, zero_below_energy=25,
             no_core_splitting=False, merge_all=None,
             gauss_width=0.0, lorentz_width=None, n_valence=None,
             core_hole_width='campbell-papp'):

        """ core energies """
        ktypes, type_map = self.core_state_types
        n_ktypes = len(ktypes)
        groups = [np.where(type_map[self.core_states_index] == i)[0].tolist() for i in range(len(ktypes))]
        if len(groups) > 1:
            if groups[0][0] > groups[1][0]:
                groups = groups[::-1]

        e_core = self.core_states_energies(no_core_splitting)
        e_core_avg = np.empty(n_ktypes)

        for i in range(n_ktypes):
            e_core_avg[i] = np.average(e_core[:,groups[i]])

        dcormin = float('inf')
        border = len(groups[0])
        if n_ktypes > 1:
            for i in e_core:
                dcormin = min(dcormin, np.min(np.abs(i[:border, np.newaxis] - i[border:])))

        """ energies """
        energies = self.energies()
        order=np.argsort(energies)
        senergies = energies[order]

        if senergies[-1] > 650:
            n_valence = 0

        """ replace energies below e-fermi and add e-fermi energy """
        fidx = np.searchsorted(senergies, 0.)
        if fidx<=0 or fidx>=len(senergies):
            raise ValueError('Fermi energy not in the range of energies')

        if abs(senergies[fidx]) < abs(senergies[fidx-1]):
            nef = fidx
        else:
            nef = fidx - 1

        if interpolate_to_fermi != False:
            if abs(senergies[nef]) < interpolation_threshold:
                interpolate_to_fermi = False

        def make_slice(order):
            """ transform order array to slice, if possible - check for subsequent numbers """
            if len(order) ==0:
                return slice(0,0)
            if np.all(order[:-1] < order[1:]) and 2*np.sum(order) == (order[-1] - order[0] + 1) * (order[0]+order[-1]):
                return slice(order[0], order[-1]+1)
            return order

        if zero_below or interpolate_to_fermi:
            if interpolate_to_fermi:
                upper_half_idx_source = upper_half_idx = nef
                #weights for interpolation to fermi
                if ((nef == 0) or (senergies[nef] <= 0.0)) and (nef < len(senergies) - 1):
                    interp_idx = nef + 1
                else:
                    interp_idx = nef - 1
                x1,x2 = senergies[nef], senergies[interp_idx]
                df = x2 - x1
                w1,w2 = -x2/df, x1/df
                my_fermi = 0.
            else:
                upper_half_idx_source = upper_half_idx = nef
                my_fermi = senergies[nef]
            if zero_below:
                senergies = np.concatenate([
                              abs(zero_below_energy) * (np.linspace(-1., 0., zero_below_num, endpoint=False)**3),
                              senergies[upper_half_idx:]
                              ])
                upper_half_idx = zero_below_num

            if interpolate_to_fermi:
                 senergies[upper_half_idx] = 0
            low_order = make_slice(order[:fidx])
            high_order = make_slice(order[upper_half_idx_source:])
        else:
            low_order = make_slice(order)
            high_order = False

        size = len(senergies)
                # Energy, Group, State, Polarization
        STATES = self.GROUPS[0].CORE_STATES
        n_states = STATES[-1]['ICST']
        n_pol = len(self.NPOL())

        shape = (size, len(self.GROUPS), n_states, n_pol)
        source_shape = (len(energies), len(self.GROUPS), n_states, n_pol)

        def read_rdt(name):
            """ Fill the below fermi with either zero or data
            fill above the fermi with data
            And finally, the middle fermi-energy point set as interpolation.
            """
            source = np.empty(source_shape)
            for i,g in enumerate(self.GROUPS.values()):
                for j,d in enumerate(g.DATA.values()):
                  source[j,i,] = d.DATA[name].real.reshape((-1,n_pol))
            source = source[order]

            out = np.empty(shape)
            if zero_below:
                out[:upper_half_idx] = 0.
                out[upper_half_idx:] = source[upper_half_idx_source:]
            else:
                out[:] = source
                if interpolate_to_fermi:
                    out[upper_half_idx] = w1 * source[nef] + w2 * source[interp_idx]
            return out

        rd = read_rdt('ABSRTA')
        rt = read_rdt('ABSRATE')

        if n_pol == 3:
            #linear dichroism - transform it to handle it as circular
            rd[:,:,:,1] = 0.5*(rd[:,:,:,0] + rd[:,:,:,1])
            rd=rd[:,:,:,1:]
            rt[:,:,:,1] = 0.5*(rt[:,:,:,0] + rt[:,:,:,1])
            rt=rt[:,:,:,1:]
            n_pol=2

        # If caller didn't provide explicit lorentz_width, compute using requested dataset
        if lorentz_width is None:
            lorentz_width = self.lorentz_width(core_hole_width)

        broaden = self.broadener(senergies, gauss_width, lorentz_width, n_valence)
        borders = [0, border, n_states]
        for i in range(n_ktypes):
            rd[:,:,borders[i]:borders[i+1]] = broaden(rd[:,:,borders[i]:borders[i+1]], i)
            rt[:,:,borders[i]:borders[i+1]] = broaden(rt[:,:,borders[i]:borders[i+1]], i)

        if merge_all is None:
            merge_all = abs(e_core_avg[0] - e_core_avg[-1]) < (senergies[-1] - senergies[0])
        if n_ktypes <= 1 or no_core_splitting:
            merge_all = False

        """
        Výpočet SD
        """

        def shift_by_e_core_avg_diff(data, energies=senergies, borders=borders, out=None, ref_avg=None, clip=False):
            """
            Interpolate `data` shifted by core-energy differences onto `energies`.

            Parameters
            - data: array shaped (len(senergies), groups, states, n_pol)
            - energies: target energy grid where the shifted spectra are evaluated
            - borders: grouping borders as used elsewhere in this function
            - out: optional pre-allocated output array; if None a new array is created
            - out: optional pre-allocated output array; if None a new array is created.

            Returns the output array (the same object passed as `out` when provided).
            """
            # determine evaluation grid and output view
            N = len(energies)
            if out is None:
                out = np.empty((N,) + data.shape[1:])

            eval_energies = energies
            out_view = out

            for i in range(e_core.shape[0]):
                for k in range(n_ktypes):
                    for j in range(borders[k], borders[k+1]):
                        if ref_avg is None:
                            de = e_core_avg[k] - e_core[i, j]
                        else:
                            de = ref_avg - e_core[i, j]
                        shifted = eval_energies - de
                        for p in range(n_pol):
                            f = interp1d(senergies, data[:, i, j, p], kind='cubic', fill_value="extrapolate")
                            out_view[:, i, j, p] = f(shifted)
                            # values corresponding to evaluation points outside source grid
                            if clip:
                                out_view[shifted < senergies[0], i, j, p] = 0.0
                                out_view[shifted > senergies[-1], i, j, p] = data[-1, i, j, p]

            return out


        if merge_all:
            ie = np.searchsorted(senergies, dcormin, side='right')
            if ie >= len(senergies):
                raise ValueError("No E(i) > Dcormin found")
            shifted_energies = np.concatenate((
                senergies[:ie],
                senergies[upper_half_idx:]-my_fermi + dcormin
            ))
            kref = np.argmax(e_core_avg)
            avg = e_core_avg[kref]

            # Use shift_by_e_core_avg_diff to fill ad/at on the merged energy grid
            # pass ref_avg so all states are shifted to the same reference (kref)
            ad = shift_by_e_core_avg_diff(rd, energies=shifted_energies, borders=borders, ref_avg=avg, clip=True)
            at = shift_by_e_core_avg_diff(rt, energies=shifted_energies, borders=borders, ref_avg=avg, clip=True)

            df = avg - e_core
            ebot = shifted_energies[0] + np.max(df)
            iebot = np.searchsorted(shifted_energies, ebot, side='right')
            etop = senergies[-1] + np.min(df)
            ietop = np.searchsorted(shifted_energies, etop, side='left')
            shifted_energies = shifted_energies[iebot:ietop]

            ad = ad[iebot:ietop]
            at = at[iebot:ietop]
            # Recompute per-ktype `sd` on the merged energy grid (preserve per-ktype shifts)
            sd = shift_by_e_core_avg_diff(rd, energies=shifted_energies, borders=borders, clip=True)
            senergies = shifted_energies
        else:
            sd = shift_by_e_core_avg_diff(rd)
            ad = sd.copy()
            at = shift_by_e_core_avg_diff(rt)



        # ------ combine spectra with common KAPPA for core state
        shape = ad.shape[:2] + (n_ktypes, n_pol)
        rd = np.empty(shape)
        rt = np.empty(shape)
        for i in range(n_ktypes):
            rd[:,:,i] = ad[:,:,borders[i]:borders[i+1]].sum(axis=2)
            rt[:,:,i] = at[:,:,borders[i]:borders[i+1]].sum(axis=2)

        index = [i.IT()-1 for i in self.GROUPS.values()]

        concealment = np.asarray(self.TYPES()['CONC']).ravel()[index]
        nqt = np.asarray(self.TYPES()['NAT']).ravel()[index]
        weight = (concealment * nqt)[None, :, None, None]

        def sum_groups(data):
            data *= weight
            return data.sum(axis=1)

        if merge_all:
            rd = rd.sum(axis = 2)[:,:,None]
            rt = rt.sum(axis = 2)[:,:,None]

        rd = sum_groups(rd)
        rt = sum_groups(rt)

        tt = 0.5*(rd[:,:,-1] + rd[:,:,0]).T
        td = 0.5*(rd[:,:,-1] - rd[:,:,0]).T

        out = {
            'POLARIZATION': tt,
            'DIFFERENCE': td,
            'ENERGY': senergies
        }

        if n_ktypes > 1:
            _sd = np.empty(shape)
            for i in range(n_ktypes):
                _sd[:,:,i]  = sd[:,:,borders[i]:borders[i+1]].sum(axis=2)
            sd = sum_groups(_sd)  #now (energy, ktype, n_pol) shape

            yd = 0.5*( sd[:,:,1] - sd[:,:,0])
            spin= 3*(yd[:,1] - 2*yd[:,0])
            orb=  2*(yd[:,1] + yd[:, 0])
            out['SPIN'] = spin
            out['ORBIT'] = orb

        self._data = out
        return out

        """
        Czech notes - mimics the XBAND plot
        //    NCXRAY(ITXG) = principal quantum number n of the absorbing core orbital.
        //    LCXRAY(ITXG) = orbital angular momentum l of that core orbital.
        TODO * CORE STATES se opakuje, ale N a L musí být stejné pro všechny skupiny
        ?? asi OK   * očísluje grupy tak, že jde číslem grupy chodit rovnou
                    do tabulky typů
        OK * vytvoří seřazenou tabulku indexů  energií E
                   - ty by měly být stejné pro všechny grupy
        OK * odečte od CORE_STATES.E(Ry) a ENERGY[0] Fermi enegy (EFERMI)
        OK * najde energii nejblíže k FE (NEF je index, DIST vzdálenost)

        OK * pokud DIST > 1E-4 tak provede lineární interpolaci do EFERMI
          z dvou nejbližších bodů (starý kód přepíše ten nejbližší bod)
    energies udělány

        OK * pro XAS/CHIXAS
              - vyhodí cokoli pod Fermi level
              - pokud ne NOTABEXT tak přidá samé nuly na n=50 energií pod Fermi level
                do souřadnic - 25Ev * (i/n)**3
        ??? * pokud je rozdíl nějakých energií menší nežli 1e-6, vynásobí tu energii 1.00001
        * rozdělí CST podle KAP, mohou být jen dva typy - pole KTYP
        OK - pokud poslední energie > 50,  nastav NVAL=0, jinak NVAL=number_of_valence_electrons_of_the_first_atom (může být předefinováno)

        OK * Načte wcorehole(atomic_number, N(1),L(1),IK=1,2 pro počet typů)

        OK * core-spliting (vypnutelné parametrem) - pokud je vypnuto
          CORE_STATES[i,j].E(Ry) =  CORE_STATES[0,0].E(Ry)

        OK * ECOREAV(KTYP)  = avg(  CORE_STATES[i,j].E(Ry) for given KTYP)
          - pokud ECOREAV(i) - ECOREAV(j) < E(1) - E(n) (v abs)
            nastav MERGE_ALL

        OK * dormin, dcormax min a max abs(rozdílů) všech CORE_STATES[i,j].E(Ry) pro core
            s ROZDILNYM KTYP. Asi jde udělat tak, že se najde min a max energií pro každý
            ktyp a pak se to odečte.

        OK  seřadí data v RD/RT (skupiny) podle energií E
          RD/RT jsou  REAL(ABSRTA), REAL(ABSRATE)
        OK broaodening RD a RT pro datasety dle energií
          íce způsobů, dopracovat
        * SD = transformace podle ECOREAV, viz line 571
              - lagrange interpolace, kdy se posunou hodnoty o ECOREAV(K) - ECORE(ITXG,ICST)
        * Pokud ne Meerge ALL, tak
          - AD, AT je stejná interpolace RD a RT, jako výš
        * Pokud Merge ALL - merging spectra, see line 620 and below

        * sum AD a AT according to the KTYP to RD, RT
          sum SD according to the KTYP (to SD)
        * 845: (*not XMO/XRS) - xmo: magneto-optics, xrs: resonant-scattering
                 over ITXG, NE
                         SUMK(K) = SUMK(K) + 0.5*(E(I)-E(I-1))
                           *(RD(I,ITXG,K,NPOL)-RD(I,ITXG,K,1))
        * Pokud MergeALL, sum RD/RT over K
        """

class RATDefinition(OutputFileDefinition):

    result_class = RATOutputFile

def G(name, **kwargs):
    def get(self, option):
        return self.option._container.data[name]

    return GeneratedValueDefinition(name, get, **kwargs)


def generate_data_calling_function(fn):

    @add_to_signature(fn)
    def wrapped(option, _inside_plot=False,
                interpolate_to_fermi=True, interpolation_threshold=1e-4,
                zero_below=True, zero_below_num=50, zero_below_energy=25,
                no_core_splitting=False, merge_all=None,
                gauss_width=0.0, lorentz_width=None, n_valence=None,
                core_hole_width='campbell-papp',
                *args, **kwargs):

        if not _inside_plot:
            option._container.generate_data(interpolate_to_fermi, interpolation_threshold,
                                            zero_below, zero_below_num, zero_below_energy,
                                            no_core_splitting, merge_all,
                                            gauss_width, lorentz_width, n_valence,
                                            core_hole_width)

        fn(option, *args, _inside_plot=_inside_plot, **kwargs)

    return wrapped


def draw_plot(axis, x, y, title, **kwargs):
      axis.set_title(title)
      axis.set_xlabel(r'$E-E_{\rm F}$ (eV)')
      axis.set_ylabel(r'XRAS')
      set_up_common_plot(axis,**kwargs)
      axis.plot(x, y)


def plot(title):

    @generate_data_calling_function
    def plot(option, axis=None, _inside_plot=False, **kwargs):
        def just_plot(plotter):
            for i, d in enumerate(data):
                if len(data)>1:
                    tit = f"{title} - type {i+1}"
                else:
                    tit = title
                with plotter(title) as axis:
                    draw_plot(axis, energy, d, tit, **kwargs)

        c = option._container
        data=c.data
        energy = data['ENERGY']
        data = data[option.name]
        if len(data.shape) == 1:
            data=data[None]

        if _inside_plot:
            just_plot(_inside_plot.new_axis)
        elif len(data.shape) and len(data) > 1:
            kwargs['layout'] = (2,1)
            with Multiplot(**kwargs) as mp:
                  just_plot(mp.new_axis)
        else:
           just_plot(single_plot)

    return plot




def create_definition():

    definition = create_output_file_definition(Keyword('RXAS'), [
      G('POLARIZATION', plot=plot('Polarization averaged spectra')),
      G('DIFFERENCE', plot=plot('Difference spectra')),
      G('SPIN', plot=plot('Spin')),
      G('ORBIT', plot=plot('Orbit')),
      G('ENERGY'),

      V('NTXRSGRP', int),
      V('IDIPOL', int),
      V('NPOL', Array(Char(), write_length=True)),
      OutputFileSectionDefinition('GROUPS', [
        V('SPECTRUM', str),
        V('IT', int),
        V('NCXRAY', int),
        V('LCXRAY', int),
        BlankSeparator(),
        BlankSeparator(),
        SeparatorDefinition(KeywordSeparator('CORE STATES :')),
        BlankSeparator(),
        V('NCST', int, written_name = 'NCST:', indent=' '),
        V('NKPCOR', Array(int), written_name = 'NKPCOR:', indent=' '),
        BlankSeparator(),
        V('CORE_STATES', Table({'ICST': int, 'N':int, 'L':int,
                         'KAP': int, 'MUE': Array(int, delimiter='/'),
                         'IKM': int, 'NORM': float,
                         'E(Ry)': float, 'E(eV)': float,
                         '<SIGMA_z>': float, 'I0': int},
                          repeat_sparse=['IKM', 'NORM']
                          ), name_in_grammar=False),
        BlankSeparator(),
        V('N_ENERGIES', int, written_name='number of energies', indent=' '),
        V('IFMT', Array(int, length=2), written_name='output format IFMT', indent=' '),
        V('VUC', float),
        BlankSeparator(),
        OutputFileSectionDefinition('DATA', [
          V('ENERGY', Sequence(
              Complex(postfix=' RYD ', delimiter='', prefix=''),
              Real(postfix=' EV (REL EF)')
          ), name_in_grammar=False),
          V('DATA',Table({'ABSRTA':complex, 'ABSRTB':complex, 'ABSRATE':complex, 'ADD':str},
                          grouping=True, numbering=True, header=False,
                          groups_as_list=False
                          ),
                  name_in_grammar= False)
        ], name_in_grammar=False, repeated_delimiter='\n', is_repeated=BaseDefinition.Repeated.LIST_SECTION)
      ], is_repeated=True, name_in_grammar=False),
    ],
    name='RAT', cls=RATDefinition, info='Result of a X-ray spectroscopy'
    )

    return definition


definition = create_definition()
