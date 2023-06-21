from ..output_files_definitions import OutputFileValueDefinition as V, create_output_file_definition, OutputFileDefinition
from typing import Optional
import numpy as np

from ..output_files import Arithmetic, CommonOutputFile
from ...common.grammar_types  import unsigned, Array, Table, RestOfTheFile, NumpyArray
from ...common.generated_configuration_definitions import NumpyViewDefinition as NV
from ...visualise.plot import PlotInfo, combined_colormap, Multiplot, PlotValue

class ARPESOutputFile(CommonOutputFile, Arithmetic):

    def plot(self, layout=(2,2), figsize=(6,4), latex=True,
             filename:Optional[str]=None, show:Optional[bool]=None, dpi=600,
             **kwargs
             ):
        mp=Multiplot(layout=layout, figsize=figsize, latex=latex)
        mp.plot(self.TOTAL, **kwargs)
        mp.plot(self.UP, **kwargs)
        mp.plot(self.DOWN, **kwargs)
        mp.plot(self.POLARIZATION, **kwargs)
        mp.finish(filename, show, dpi)
        breakpoint()

    _arithmetic_values = [('RAW_DATA', (slice(None), slice(2,6)))]

    def _check_arithmetic(self, other):
        try:
            assert np.allclose(other.ENERGY(), self.ENERGY())
            assert np.allclose(other.THETA(), self.THETA())
        except AssertionError:
            raise ValueError("The outputs are not compatibile to summed/subtracted.");


class ARPESDefinition(OutputFileDefinition):

    result_class = ARPESOutputFile


def create_definition():

    pi = PlotInfo(
        axes = (('K', r'$k_{\parallel} $(${\rm \AA}^{-1}$)'),
                ('ENERGY', r'$E-E_{\rm F}$ (eV)')),
        #axes = ('K', 'ENERGY' ),
        show_zero_line = True,
        mode = 'from_zero',
        norm='log',
        vmax = PlotValue( lambda o: o._container.TOTAL().max() )
    )
    def i(j):
        return slice(None),j

    definition = create_output_file_definition('ARPES', [
      V('NT', int),
      V('NP', int),
      V('RAW_DATA', NumpyArray(written_shape=(-1,8)), name_in_grammar=False),

      NV('THETA', 'RAW_DATA', i(0), ('NE', 'NT')),
      NV('ENERGY', 'RAW_DATA', i(1), ('NE', 'NT')),
      NV('TOTAL', 'RAW_DATA', i(2), ('NE', 'NT'), info='Total intensity', plot=pi),
      NV('UP', 'RAW_DATA', i(3), ('NE', 'NT'), info='Spin up', plot=pi),
      NV('DOWN', 'RAW_DATA', i(4), ('NE', 'NT'), info='Spin down', plot=pi),
      NV('POLARIZATION', 'RAW_DATA', i(5), ('NE', 'NT'), info='Spin polarization',
         plot=pi + {
           'mode' : 'zero_centered',
           'norm' : 'lin',
           'vmax' : None
         },
        ),
      NV('K', 'RAW_DATA', i(6), ('NE', 'NT'), info='K_parallel (pi/A)'),
      NV('DETERMINANT', 'RAW_DATA', i(7), ('NE', 'NT')),
    ], cls=ARPESDefinition)
    return definition

definition = create_definition()
