from ..output_files_definitions import OutputFileValueDefinition as V, create_output_file_definition, OutputFileDefinition
from typing import Optional
import numpy as np

from ..output_files import Arithmetic, CommonOutputFile
from ...common.grammar_types import NumpyArray, Prefixed
from ...common.generated_configuration_definitions import NumpyViewDefinition as NV
from ...visualise.plot import change_default_kwargs, colormesh, Multiplot
import matplotlib.pyplot as plt


class ARPESOutputFile(CommonOutputFile, Arithmetic):

    def plot(self, layout=(2,2), figsize=(10,6), latex=True,
             filename:Optional[str]=None, show:Optional[bool]=None, dpi=600,
             **kwargs
             ):
        mp=Multiplot(layout=layout, figsize=figsize, latex=latex)
        plt.subplots_adjust(left=0.12,right=0.95,bottom=0.17,top=0.90, hspace=0.75, wspace=0.5)
        mp.plot(self.TOTAL, **kwargs)
        mp.plot(self.UP, **kwargs)
        mp.plot(self.DOWN, **kwargs)
        mp.plot(self.POLARIZATION, **kwargs)
        mp.finish(filename, show, dpi)

    _arithmetic_values = [('RAW_DATA', (slice(None), slice(2,6)))]

    def _assert_arithmetic(self, other):
         """ Check, that the file can be summed/subtracked from an other file """
         assert np.allclose(other.ENERGY(), self.ENERGY())
         assert np.allclose(other.THETA(), self.THETA())


class ARPESDefinition(OutputFileDefinition):

    result_class = ARPESOutputFile


def create_definition():

    def plot(option, **kwargs):
        c = option._container
        kw = {
            'show_zero_line' : True,
            'mode' : 'from_zero',
            'norm' : 'log',
            'vmax' : c.TOTAL().max(),
            'xlabel' : r'$k_{\parallel} $(${\rm \AA}^{-1}$)',
            'ylabel' : r'$E-E_{\rm F}$ (eV)',
        }
        kw.update(kwargs)
        colormesh(c.K(), c.ENERGY(), option(), **kw)

    def i(j):
        return slice(None),j

    definition = create_output_file_definition('ARPES', [
      V('NT', int),
      V('NP', int),
      V('COMMENT', Prefixed('#'), name_in_grammar=False),
      V('RAW_DATA', NumpyArray(written_shape=(-1,8)), name_in_grammar=False),

      NV('THETA', 'RAW_DATA', i(0), ('NE', 'NT')),
      NV('ENERGY', 'RAW_DATA', i(1), ('NE', 'NT')),
      NV('TOTAL', 'RAW_DATA', i(2), ('NE', 'NT'), info='Total intensity', plot=plot),
      NV('UP', 'RAW_DATA', i(3), ('NE', 'NT'), info='Spin up', plot=plot),
      NV('DOWN', 'RAW_DATA', i(4), ('NE', 'NT'), info='Spin down', plot=plot),
      NV('POLARIZATION', 'RAW_DATA', i(5), ('NE', 'NT'), info='Spin polarization',
         plot=change_default_kwargs(plot, mode = 'zero_centered', norm = 'lin', vmax = None)
        ),
      NV('K', 'RAW_DATA', i(6), ('NE', 'NT'), info='K_parallel (pi/A)'),
      NV('DETERMINANT', 'RAW_DATA', i(7), ('NE', 'NT')),
    ], cls=ARPESDefinition, info='ARPES (Angle-resolved photoemission spectroscopy) output.')
    return definition


definition = create_definition()
