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


    reorder = (1, 0, 2)

    def i(type):
        def index(data, c):
          """
          Returns data in the shape
          ('NE/NK1','NQ_EFF', 'NK2')
          for a given type.
          Reorder parameter then change the order of the axes to
          ('NQ_EFF', 'NE/NK1', 'NK2')

          Data structure:
            IX,Y,Z (BSF-SPOL)
                k-e: NE, type, NQ, NK  types(I,x,y,z)
                k-k: type, K1, NQ, K2, types(I,x,y,z)

           Iup,dn (BSF) two sets
                k-e: NE, type, NQ, NK   types(u,d)
                     NE, type, NQ, NK   types(I)
                k-k: NK, type, NQ, NK   types(u,d)
                     NK, NQ, NK   types(I)
          """
          nq=c.NQ_EFF()
          nk2=c.NK2()
          if c.KEYWORD() == 'BSF':
             nk1 = c.NE() if c.MODE() == 'EK-REL' else c.NK1()
             limit = nq*nk2*2
             if type:
                return data[:limit].reshape(nk1, 2, nq, nk2)[:,type]
             else:
                return data[limit:].reshape(nk1,nq,nk2)
          else:
             if c.MODE() == 'EK-REL':
                return data.reshape(c.NE(), 4, nq, nk2)[:,type]
             else:
                return data.reshape(4, c.NK1(), nq, nk2)[:,type]
        return index
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
            NV('I_UP', 'RAW_DATA', i(c,1), reorder=reorder ),
            NV('I_DOWN', 'RAW_DATA', i(c,2), reorder=reorder ),
        ],
        'BSF-SPOL': [
            NV('I_X', 'RAW_DATA', i(c,1), reorder=reorder ),
            NV('I_Y', 'RAW_DATA', i(c,2), reorder=reorder ),
            NV('I_Z', 'RAW_DATA', i(c,3), reorder=reorder ),
        ]
      }),
      NV('I', 'RAW_DATA', i(c, 0), reorder=reorder ),

    ], cls=BSFDefinition, name='BSF')

    #definition['INDKDIR'].add_grammar_hook(lambda x: x.addParseAction(lambda x: breakpoint() or x))
    return definition

definition = create_definition()
