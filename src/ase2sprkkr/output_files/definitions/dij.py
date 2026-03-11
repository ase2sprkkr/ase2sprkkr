from ..output_files_definitions import OutputFileValueDefinition as V, create_output_file_definition, OutputFileDefinition
from typing import Optional
import numpy as np
import re

from ..output_files import Arithmetic, CommonOutputFile
from ...common.grammar_types import RawData, RestOfTheFile, Table, Integer, Real, \
                                    Array, Sequence, String, LineString, NumpyArray
from ...common.generated_configuration_definitions import NumpyViewDefinition as NV, \
                                                          GeneratedValueDefinition as GV
from ...common.configuration_definitions import SeparatorDefinition
#from ...gui.plot import change_default_kwargs, colormesh, Multiplot
import matplotlib.pyplot as plt

class JXCOutputFile(CommonOutputFile):
    pass


class JXCOutputFileDefinition(OutputFileDefinition):
    result_class = JXCOutputFile

import pyparsing as pp
pp.ParserElement.verbose_stacktrace = True

def create_definition():

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
      V('TABLE_HEADER', LineString(), default_value='IT   IQ   JT    JQ   N1 N2 N3    DRX    DRY    DRZ     DR       DX_ij [meV]    DY_ij [meV]    DZ_ij [meV]', is_hidden=True, name_in_grammar=False),
      V('DATA', NumpyArray(dtype=[
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
            ('DX', float),
            ('DY', float),
            ('DZ', float),
        ]), name_in_grammar=False)
    ], info='Dzyaloshinski-Moriya couplings Dij')

    definition.__dict__['extension'] = 'dat'

    return definition


definition = create_definition()


"""
      V('NP', int),
      V('COMMENT', Prefixed('#'), name_in_grammar=False),
      V('RAW_DATA', NumpyArray(written_shape=(-1,8)), name_in_grammar=False),
      GV('MODE', lambda c,k=None: 'energy' if c.NE() > 1 else 'kx_ky'),
      *switch('MODE', {
        'energy': [
          NV('THETA', 'RAW_DATA', i(0), ('NE', 'NT')),
          NV('ENERGY', 'RAW_DATA', i(1), ('NE', 'NT')),
          NV('K', 'RAW_DATA', i(6), ('NE', 'NT'), info='K_parallel (pi/A)'),
        ],
        'kx_ky': [
          NV('KX', 'RAW_DATA', i(0), ('NP', 'NT'), info='k_x coordinate (1/Å)'),
          NV('KY', 'RAW_DATA', i(1), ('NP', 'NT'), info='k_y coordinate (1/Å)'),
        ]
      }),
      NV('TOTAL', 'RAW_DATA', ii(2), info='Total intensity', plot=change_default_kwargs(plot,title=r'Total intensity',colormap='gray')),
      NV('UP', 'RAW_DATA', ii(3), info='Spin up', plot=change_default_kwargs(plot,title=r'Spin up',colormap='gray')),
      NV('DOWN', 'RAW_DATA', ii(4), info='Spin down', plot=change_default_kwargs(plot, title=r'Spin down', colormap='gray')),
      NV('POLARIZATION', 'RAW_DATA', ii(5), info='Spin polarization',
         plot=change_default_kwargs(plot, colormap='bwr', title=r'Spin polarization')
        ),
      NV('DETERMINANT', 'RAW_DATA', determinant),
"""
