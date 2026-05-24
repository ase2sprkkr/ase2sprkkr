from ..output_files_definitions import OutputFileValueDefinition as V, OutputFileDefinition
from typing import Optional
import numpy as np
import re
import itertools
from enum import Enum
from matplotlib import cm
import warnings
from collections import namedtuple

from ..output_files import OutputFile
from ...common.decorators import cached_property, cached_class_property
from ...common.grammar_types import RawData, Table, Integer, Real, \
                                    Array, Sequence, String, LineString, NumpyArray
from ...common.configuration_definitions import SeparatorDefinition
from ...gui.plot import Multiplot

class Coordinates(Enum):
    cartesian = "cartesian"
    lattice = "lattice"


Selector = namedtuple("Selector", ["iq", "it"])

class JXCOutputFile(OutputFile):

    plot_parameters = {'exchange_radius', 'iq', 'it', 'exclude_it', 'exclude_vc', 'font_size', 'axis', 'separate_plots'}

    def is_Jij(self):
        return len(self.DATA().dtype.names) == 15

    def is_Dij(self):
        return len(self.DATA().dtype.names) == 14

    @cached_property
    def iqs(self):
        return {
            i[0]: { d[0]: (d[1], d[2]) for d in i[2] }  for i in self.OCCUPATION()
        }

    @cached_property
    def it_labels(self):
        def generator():
            for i in self.iqs.values():
                for it, (label, occ) in  i.items():
                    yield it, label

        return dict(generator())

    def iq_to_its(self, iq):
        return self.iqs[iq].keys()

    def it_to_iqs(self, it):
        return { iq for iq, types in self.iqs.items() if it in types }

    @cached_property
    def labels_to_it(self):
        return {label: it for it, label in self.it_labels.items()}

    def label_to_it(self, label):
        try:
            return self.labels_to_it[label]
        except KeyError as exc:
            raise ValueError(f"Unknown type label '{label}'. Valid labels are: {', '.join(self.labels_to_it.keys())}") from exc

    def element_to_its(self, label:str):
        return [ it for it, l in self.it_labels.items() if l.startswith(label) and (l == label or l[len(label)] == '_') ]

    @cached_class_property
    def _it_selector_regex():
        return re.compile(r'\s*([A-Z][a-z]*(_(\d+)))?\s*')

    def create_selector(self, selector=None, iq=None, it=None, exclude_it=None, exclude_vc=True):
        if selector is not None:
            return selector

        def parse_list(values):
            if isinstance(values, str):
                return values.split(',')
            if isinstance(values, int):
                return [values]
            return values

        def parse_it(values):
            if values is None or values == "":
                return ...
            values = parse_list(values)

            out = set()

            for v in values:
                if isinstance(v, int):
                    out.add(v)
                    continue
                match = self._it_selector_regex.fullmatch(v)
                if not match:
                    try:
                        index = int(v)
                        if index not in self.it_labels:
                            raise ValueError("Selector index out of range: {v}")
                        out.add(index)
                    except ValueError:
                        raise ValueError(f"Invalid selector: {v}")
                elif match.group(2):
                    out.add(self.label_to_it(match.group(1)))
                else:
                    its = self.element_to_its(match.group(1))
                    if not its:
                        raise ValueError(f'Selector {v} does not match any atomic type')
                    out.update(its)
            return out

        def parse_iq(values):
            if values is None or values == "":
                return ...
            values = parse_list(values)
            try:
                return [ int(v) for v in values ]
            except:
                raise ValueError(f"Invalid IQ selector: {values}. Must be (a comma-separated) list of integers.")

        iq = parse_iq(iq)
        it = parse_it(it)
        exclude_it = parse_it(exclude_it)

        if exclude_vc:
            exclude_vc = set(self.element_to_its('Vc'))
            if exclude_it is ...:
                exclude_it = exclude_vc
            else:
                exclude_it.update(exclude_vc)

        if exclude_it is not ...:
            exclude_it = parse_it(exclude_it)
            if it is ...:
                it = set(range(1, self.NT() + 1))
            it = it - exclude_it
        return Selector(iq, it)

    def it_selector(self, selector=None, iq=None, it=None, exclude_it=None, exclude_vc=True):
        if selector is None:
            selector = self.create_selector(iq=iq, it=it, exclude_vc=exclude_vc, exclude_it=exclude_it)
        iq, it = selector
        out = ...
        if iq is not ...:
            out = itertools.chain.from_iterable(self.iq_to_its(i) for i in iq)
        if it is not ...:
            out = it if out is ... else set(out).intersection(it)
        return out

    def iq_selector(self, selector=None, iq=None, it=None, exclude_it=None, exclude_vc=True):
        if selector is None:
            selector = self.create_selector(iq=iq, it=it, exclude_it=exclude_it, exclude_vc=exclude_vc)
        iq, it = selector

        out = self.it_selector(selector=selector)
        if out is ...:
            return ...
        return set().union(* (self.it_to_iqs(it) for it in out) )

    def filtered_data(self, selector=None, iq=None, it=None, exclude_it=None, exclude_vc=True, exchange_radius=4.0):
        its = self.it_selector(selector=selector, iq=iq, it=it, exclude_it=exclude_it, exclude_vc=exclude_vc)
        data = self.DATA()

        if its is ...:
            mask = None
        else:
            its = list(its)
            mask = np.isin(data["IT"], its)
            mask &= np.isin(data["JT"], its)

        if exchange_radius is not None:
            if mask is not None:
                mask &= data['DR'] <= exchange_radius
            else:
                mask = data['DR'] <= exchange_radius

        if mask is not None:
            data = data[mask]
        return data

    def _spin_moments(self):
        atomic_types = set()
        moments = []
        for site in self.atoms.sites:
            for atomic_type in site.occupation:
                if atomic_type in atomic_types:
                    continue
                atomic_types.add(atomic_type)
                moment = getattr(atomic_type.moments, 'spin_moment', None)
                if moment is None:
                    return None
                moments.append(moment)
        return moments

    @staticmethod
    def _dij_components(component_axis):
        components = {
            'x': (['DX'], ['x']),
            'y': (['DY'], ['y']),
            'z': (['DZ'], ['z']),
            'all': (['DX', 'DY', 'DZ'], ['x', 'y', 'z']),
        }
        component_axis = { 0: 'x', 1: 'y', 2: 'z', None: 'all' }.get(component_axis, component_axis)
        try:
            return components[component_axis]
        except KeyError as exc:
            raise ValueError(f"Invalid DMI axis '{component_axis}'. Use one of: all, x, y, z.") from exc

    def plot(self, layout=2, figsize=None, latex=None,
             filename:Optional[str]=None, show:Optional[bool]=None, dpi=300, label_spacing=0.3,
             selector=None,
             iq=None, it=None, exclude_it=None, exclude_vc=True, exchange_radius=4.0,
             font_size=10, axis='all',
             separate_plots=False, layout_kind = 'constrained',
             **kwargs):

        def _resolve_layout(count, values):
            if isinstance(layout, int):
                if values > 1:
                    return (count, values)
                cols = min(layout, values)
                rows = count // cols + (1 if count % cols else 0)
                return (max(rows, 1), cols)
            return layout

        labels = self.it_labels

        data = self.filtered_data(selector=selector, iq=iq, it=it, exclude_it=exclude_it, exclude_vc=exclude_vc, exchange_radius=exchange_radius)
        type_indexes = sorted(set(labels.keys()))
        if not type_indexes or not len(data):
            return False
        is_jij = self.is_Jij()
        x = data['DR']

        if is_jij:
            label = lambda partner_index: f'${labels[type_index]}$-${labels[partner_index]}$'
            value_label = lambda: r'$J_{ij}$'
            name = lambda: f'Jij_{labels[type_index]}'
            y = data[['JXX']]

            spin_mom = self._spin_moments()
            if spin_mom is None:
                warnings.warn("Spin moments not found.")
            elif len(spin_mom) != self.NT():
                warnings.warn("Spin moments count does not match number of types - the potential do not matches the Jij data.")
                spin_mom = None

        else:
            index, axis_labels = self._dij_components(axis)
            label = lambda partner_index: f'{labels[type_index]}-{labels[partner_index]}' #$_{axis_labels[colindex]}$'
            value_label = lambda: rf'$D_{{ij}}^{{{axis_labels[colindex]}}}$'
            name = lambda: f'Dij_{labels[type_index]}'
            y = data[index]
            spin_mom = None

        extremum = max( ( np.abs(y[col]).max() for col in y.dtype.names ) )
        extremum += max(0.1, extremum * 0.1)
        colorspace = cm.Paired(np.linspace(0, 1, 2 * self.NT() + 2))

        layout = _resolve_layout(len(type_indexes), len(y.dtype.names))

        def plot(axis):
            colors = iter(colorspace)
            for partner_index in type_indexes:
                mask = other == partner_index
                if not np.any(mask):
                        continue
                xxx = xx[mask]
                yyy = col[mask]
                if spin_mom is not None:
                    sign_ij = np.sign(spin_mom[type_index-1]) * np.sign(spin_mom[partner_index-1])
                    yyy*=sign_ij
                #cycle over DX,DY,DZ (or just plot one line for all other cases)
                color = next(colors)
                axis.plot(xxx, yyy, alpha=0.75, lw=2.0, color=color, label=label(partner_index))
                axis.scatter(xxx, yyy, color=color, alpha=0.75, s=45, lw=1.0, edgecolor='black')
            axis.set_xlabel(r'$r_{ij}/a_{\mathrm{lat}}$', fontsize=font_size)
            axis.set_ylabel(value_label() + r' $[meV]$', fontsize=font_size)
            axis.tick_params(axis='x', colors='black', labelsize=font_size)
            axis.tick_params(axis='y', colors='black', labelsize=font_size)
            axis.axhline(0, color='black', linestyle='--')
            if label_spacing:
                axis.legend(fontsize=font_size - 2, loc='best', labelspacing=label_spacing)
            axis.grid(False)
            axis.set_ylim(-extremum, extremum)

        with Multiplot(layout=layout, figsize=figsize, latex=latex,
                           filename=filename, show=show, dpi=dpi,
                           separate_plots=separate_plots, layout_kind=layout_kind,
                           **kwargs) as mp:

                for type_index in type_indexes:
                    mask = data['IT'] == type_index
                    if not np.any(mask):
                        continue
                    xx = x[mask]
                    yy = y[mask]
                    other = data['JT'][mask]

                    for colindex, cname in enumerate(y.dtype.names):
                        col = yy[cname]

                        mp.plot(
                            self,
                            name=name(),
                            plot_function=plot,
                        )
        return True

    def write_uppasd_file(self, file_name=None, directory=None,
                          selector=None,
                          it=None, iq=None, exclude_it=None, exclude_vc=True, exchange_radius=None,
                          coordinates: Coordinates = Coordinates.lattice):
        """ Write Jij or Dij data to file in format suitable for UppAsd programm.
        If the filename is not given, the standard name for the given file type
        will be used.
        """

        if self.is_Jij():
            from ...bindings.uppasd import write_jfile as write
        else:
            from ...bindings.uppasd import write_dmfile as write
        write(self, file_name, directory=directory,
            selector=selector,
            iq=iq, it=it, exclude_it=exclude_it, exclude_vc=exclude_vc,  exchange_radius=exchange_radius, coordinates=coordinates)



    @classmethod
    def from_atoms(cls, atoms, data = None):
        out = cls(definition=definition)

        its = {}
        def it(atomic_type):
            if atomic_type in its:
                return its[atomic_type]
            index = len(its) + 1
            its[atomic_type] = index
            return index

        first = {}
        duplicates = {}

        for site in atoms.sites:
            for atomic_type in site.occupation:
                if atomic_type.symbol in first:
                    if not atomic_type.symbol in duplicates:
                        duplicates[atomic_type.symbol] =  { first[atomic_type.symbol] : 1}
                    duplicates[atomic_type.symbol][atomic_type] = len(duplicates[atomic_type.symbol]) + 1
                else:
                    first[atomic_type.symbol] = atomic_type

        def label(atomic_type):
            dups = duplicates.get(atomic_type.symbol, None)
            if dups is not None:
                return f'{atomic_type.symbol}_{dups[atomic_type]}'
            return atomic_type.symbol

        out.OCCUPATION = [
            ( i + 1, len(site.occupation),  [  ( it(atomic_type), label(atomic_type), occ) for atomic_type, occ in site.occupation.items() ] )
            for i, site in enumerate(atoms.sites)
        ]

        if data:
            out.DATA = data

        return out



class JXCOutputFileDefinition(OutputFileDefinition):
    result_class = JXCOutputFile

import pyparsing as pp
pp.ParserElement.verbose_stacktrace = True

def create_definition():

    def table_header(c):
        if c.is_Jij():
          return 'IT   IQ   JT    JQ   N1 N2 N3    DRX    DRY    DRZ     DR      J_xx [meV]     J_yy [meV]      J_xy [meV]     J_yx [meV]'
        else:
          return 'IT   IQ   JT    JQ   N1 N2 N3    DRX    DRY    DRZ     DR       DX_ij [meV]    DY_ij [meV]    DZ_ij [meV]'

    shared_dtype = [
            ('IT', int),
            ('IQ', int),
            ('JT', int),
            ('JQ', int),
            ('N1', int),
            ('N2', int),
            ('N3', int),
            ('DRX', float),
            ('DRY', float),
            ('DRZ', float),
            ('DR', float),
    ]

    definition = JXCOutputFileDefinition('JXC', [
      V('HEADER', RawData(ends_with=re.compile('\n[ \t]*number of sites'), default_value="""

 *******************************************************************************
                                   <XCPLTENSOR>:
                      Dzyaloshinski-Moriya couplings  Dij
                      according to Phys. Rev. B 79, 045209 (2009)
 *******************************************************************************

"""), name_in_grammar=False),
      V('NQ', int, written_name="number of sites   NQ", delimiter=' = ', delimiter_grammar='='),#, indent=10*" "),
      V('NT', int, written_name="number of types   NT", delimiter=' = ', delimiter_grammar='=', indent=10*" "),
      SeparatorDefinition('                              site occupation:'),
      V('OCCUPATION', Table(IQ=Integer(prefix="IQ = ", format="{:>3}"), NOQ=int,
                            DATA=Array(Sequence(int,
                                       String(prefix='-', format='{:>10}'),
                                       Real(prefix='x = ', format='{:6.3f}')),
                                       prefix='IT:'),
                            header=False
                            )),
      V('TABLE_HEADER', LineString(), default_value_from_container=table_header,
        is_hidden=True, name_in_grammar=False),
      V('DATA', NumpyArray(dtypes=[
            shared_dtype + [
                ('JXX', float),
                ('JYY', float),
                ('JXY', float),
                ('JYX', float),
            ],
            shared_dtype + [
                ('DX', float),
                ('DY', float),
                ('DZ', float),
            ],
        ]), name_in_grammar=False)
    ], info='Dzyaloshinski-Moriya couplings Dij or Jij')

    definition.__dict__['extension'] = 'dat'

    return definition


definition = create_definition()
