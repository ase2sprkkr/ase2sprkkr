from ..output_files_definitions import OutputFileValueDefinition as V, create_output_file_definition, OutputFileDefinition,Separator
from typing import Optional
import numpy as np

from ..output_files import Arithmetic, CommonOutputFile
from ...common.grammar_types  import unsigned, Array, Table, RestOfTheFile, NumpyArray, Keyword
from ...common.generated_configuration_definitions import NumpyViewDefinition as NV, GeneratedValueDefinition as GV
from ...common.configuration_definitions import gather,switch
from ...visualise.plot import PlotInfo, combined_colormap, Multiplot, PlotValue
import matplotlib.pyplot as plt

class BSFOutputFile(CommonOutputFile, Arithmetic):

    """
    Output file for Bloch spectral functions
    """

class BSFDefinition(OutputFileDefinition):

    result_class = BSFOutputFile


def create_definition():

    def i(c, i):
        if i == -1:
           i=2 if  c.MODE() == 'BSF' else 0
        ln = c.NE()*c.NQ_EFF()*c.NK()
        start=ln*i
        return slice(start,start+ln)

    reorder = (1, 0, 2)

    definition = create_output_file_definition(Keyword('BSF', 'BSF-SPOL'), [
      Separator(),
      V('DATASET', str, written_name='#DATASET'),
      V('MODE', Keyword('EK-REL', 'CONST-E')),
      *switch('MODE', {
          'EK-REL' : [
              V('NE2', int, written_name='NE', info="Number of energies (the second axis)"),
              V('NK', int, info="Number of K points (the same as NK2)"),
              Separator(),
              *gather(
                V('EMIN', float),
                V('EMAX', float)
              ),
              Separator(),
              V('NKDIR', int),
              V('LBLKDIR', NumpyArray(written_shape=(-1,1), shape=(-1,), lines='NKDIR', dtype='line'),
                           name_in_grammar=False),
              Separator(),
              V('INDKDIR', int, is_repeated=True),
              Separator(),
              V('NK2', int, written_name='NK', info="Number of K points (the last axis)"),

              GV('SHAPE', lambda c: c['NE']()),
              V('K', NumpyArray(written_shape=(-1,1), shape=(-1,), lines='NK'), name_in_grammar=False),
              Separator(),
              V('E', Array(float), init_by_default=True, default_value=lambda o:
                 np.linspace(o._container.EMIN(),o._container.EMAX(),data.NE), is_stored=False),
            ],
         'CONST-E' : [
              V('NK1', int, info='Number of K points (the first axis)'),
              V('NK2', int, info='Number of K points (the first axis)'),
              V('ERYD', Array(float, length=2)),
              V('NK1_1', int, written_name='NK1', is_hidden=True),
              V('VECK1', Array(float, length=3)),
              V('NK2_1', int, written_name='NK2', is_hidden=True),
              V('VECK2', Array(float, length=3)),
              Separator(),
              V('K1', Array(float), init_by_default=True, default_value=lambda o:
                 np.linspace(0,np.linalg.norm(data.VECK1),data.NK1), is_stored=False),
              V('K2', Array(float), init_by_default=True, default_value=lambda o:
                 np.linspace(0,np.linalg.norm(data.VECK2),data.NK2), is_stored=False),
          ]}
      ),
      V('RAW_DATA', NumpyArray(written_shape=(-1,1), shape=(-1,)), name_in_grammar=False),
      *switch('KEYWORD', {
        'BSF' : [
            NV('I_UP', 'RAW_DATA', lambda c: i(c,0), ('SHAPE','NQ_EFF', 'NK2'), reorder=reorder ),
            NV('I_DOWN', 'RAW_DATA', lambda c: i(c,1), ('SHAPE','NQ_EFF', 'NK2'), reorder=reorder ),
        ],
        'BSF-SPOL': [
            NV('I_X', 'RAW_DATA', lambda c: i(c,1), ('SHAPE','NQ_EFF', 'NK2'), reorder=reorder ),
            NV('I_Y', 'RAW_DATA', lambda c: i(c,2), ('SHAPE','NQ_EFF', 'NK2'), reorder=reorder ),
            NV('I_Z', 'RAW_DATA', lambda c: i(c,3), ('SHAPE','NQ_EFF', 'NK2'), reorder=reorder ),
        ]
      }),
      NV('I', 'RAW_DATA', lambda c: i(c, -1), ('SHAPE','NQ_EFF', 'NK2'), reorder=reorder ),

    ], cls=BSFDefinition, name='BSF')

    #definition['INDKDIR'].add_grammar_hook(lambda x: x.addParseAction(lambda x: breakpoint() or x))
    return definition

definition = create_definition()
