# coding: utf-8
# vim: set et ts=4 sw=4 nu fdm=indent cc=80:

"""
Module sprkkr_inputs
===============

This function contains input templates for various input_parameters
New template include into dict_inp_data.
It returns list of strings which can be proceeded by
load_parameters
Please be carefull, each keyword in given section has to be
on the new line
.. code-block:: python

    from sprkkr.data import sprkkr_inputs
    task='scf'
    print(get_sprkkr_input(task))

"""


dict_inp_data={}
########################################################
# INPUT TEMPLATE SCF
dict_inp_data['SCF']='''
CONTROL DATASET=case
        ADSI=SCF
        POTFIL=case.pot
        KRWS=1

TAU BZINT=POINTS
    NKTAB=250

ENERGY  GRID={5}
      NE={32}
      ImE=0.00
      EMIN=-0.2

SCF NITER=200
    MIX=0.20000
    VXC=VWN
    EFGUESS=0.7
    TOL=0.000010
    ISTBRY=1


SITES NL={3}
'''
########################################################
# INPUT TEMPLATE DOS
dict_inp_data['DOS']='''
CONTROL DATASET=case
        ADSI=DOS
        POTFIL=case.pot_new


TAU BZINT=POINTS
    NKTAB=250

ENERGY  GRID={3}
      NE={300}
      ImE=0.01
      EMIN=-0.2
      EMAX= 1.0

TASK DOS

SITES NL={3}
'''

########################################################
# INPUT TEMPLATE ARPES
dict_inp_data['ARPES']='''
CONTROL  DATASET     = case
         ADSI        = ARPES
         POTFIL      = case.pot_new
         PRINT = 0

TAU      BZINT= POINTS
         NKTAB= 250

ENERGY   GRID={1}
         NE={300}
      ImE=0.0 Ry
      EMINEV=-10.
      EMAXEV=10.0
      EWORK_EV=4.2
      IMV_INI_EV=0.05
      IMV_FIN_EV=5.

SITES 	 NL=4

TASK     ARPES
         STRVER=0
         IQ_AT_SURF=2
         MILLER_HKL={0,0,1}
         NTMP=1
         TMPMIN=11.
         CTMPMAX=11.
         VIBRA CNVIBRA=14

SPEC_PH THETA=45.
         PHI=0.
         POL_P=P
         EPHOT= 6675.

SPEC_STR N_LAYDBL={10,10}
         NLAT_G_VEC=57
         N_LAYER=120
         SURF_BAR={0.25,0.25} #TRANSP_BAR

SPEC_EL THETA={-10.,10.}
         PHI={0.,0.}
         NT=300
         NP=1
         POL_E=PZ
'''
########################################################
# INPUT TEMPLATE SCF
dict_inp_data['PHAGEN']='''
CONTROL DATASET=case
        ADSI=PHAGEN
        POTFIL=case.pot_new

TAU BZINT=POINTS
    NKTAB=250

ENERGY  GRID={5}
      NE={32}
      ImE=0.00
      EMIN=-0.2

SCF NITER=200
    MIX=0.20000
    VXC=VWN
    EFGUESS=0.7
    TOL=0.000010
    ISTBRY=1

TASK PHAGEN
'''
########################################################
########################################################
def str2lines(string):
    li = list(string.split("\n"))
    return li

def get_sprkkr_input(task):

    lines='None'

    if input_parameters.upper() == 'ARPES':
        lines=dict_inp_data['ARPES']
    elif input_parameters.upper() == 'SCF':
        lines=dict_inp_data['SCF']
    elif input_parameters.upper() == 'DOS':
        lines=dict_inp_data['DOS']
    elif input_parameters.upper() == 'PHAGEN':
        lines=dict_inp_data['PHAGEN']
    else:
        lines='None'
    return  str2lines(lines)
