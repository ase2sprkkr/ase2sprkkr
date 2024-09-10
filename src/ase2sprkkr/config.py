""" This module contains configuration, that could be changed, preferrably
by .config/ase2sprkkr/__init__.py file"""

import os
from .common.grammar_types import CustomMixed, QString, Array, Bool, Keyword, Integer
from .common.container_definitions import SectionDefinition
from .common.value_definitions import ValueDefinition as V
import functools
import warnings
import shutil


class Section(SectionDefinition):
    info_in_data_description = True


def _get_suffix(*_):
    return os.environ.get('SPRKKR_EXECUTABLE_SUFFIX','')


@functools.lru_cache
def find_default_mpi_runner():
   for r in [ 'mpirun', 'mpirun.opmpirun', 'mpirun.mpich' ]:
       if shutil.which(r):
           return [ r ]
   return False


@functools.lru_cache
def get_default_mpi_runner():

   out = find_default_mpi_runner()
   if out:
       return out
   if config.mpi_warning():
       warnings.warn("No MPI runner found. Disabling MPI!!!")


def mpi_runner(mpi):
    """ Return a shell command to execute a mpi task.

    Parameters
    ----------
    mpi_runner: Union[bool,str,list,int]


      - If True is given, return the default mpi-runner
      - If False is given, no mpi-runner is returned.
      - If 'auto' is given, it is the same as True, however
           * no warning is given if no mpi is found
           * MPI is not used, if only one CPU is available
      - If a string is given, it is interpreted as a list of one item
      - If a list (of strings) is given, the user specified its own runner, use it as is
        as the parameters for subprocess.run.
      - If an integer is given, it is interpreted as the number of
        processes: the default mpi-runner is used, and the parameters
        to specify the number of processes.

    Return
    ------
    mpi_runner: list
      List of strings with the executable and its parameters, e.g.

      ::

          ['mpirun', '-np', '4']
    """
    if mpi is None:
       mpi=config.running.mpi()
    if mpi is False:
       return None
    if mpi is True:
       mpi=find_default_mpi_runner()
    if isinstance(mpi, list):
        return mpi
    if isinstance(mpi, str):
        if mpi == 'auto':
            if hasattr(os, 'sched_getaffinity') and len(os.sched_getaffinity(0))==1:
                return None
            return find_default_mpi_runner()
        return [ mpi ]
    if isinstance(mpi, int):
       return find_default_mpi_runner() + ['-np', str(mpi)]
    return mpi


""" The definition of ASE2SPRKKR configuration """
definition = Section('config', [

  Section('running', [
    V('empty_spheres', CustomMixed(Bool.I, Keyword('auto')), default_value='auto', info="Run empty spheres finding before calculation? Default value ``auto`` means only for SCF calculations not containing any vacuum atom."),
    V('print_output', CustomMixed(Bool.I, Keyword('info')), default_value='info', info="Print output of SPRKKR calculation. Default value ``info`` prints only short info each iteration."),
    V('mpi', CustomMixed(Bool, Array(QString.I), Integer.I), is_optional=True, default_value=None,
             info='Use mpi for calculation? List of strings means yes, use the given strings as mpi runner and its params (e.g. [ "mpirun", "-n", "4" ]). Default None means try to autodetect. Integer number means use the standard runner with a given number of processes.'),
    V('mpi_warning', True, info='Warn, if no MPI is found.')
  ], info='Default values for SPRKKR calculator parameters.'),

  Section('executables', [
    V('suffix', QString.I,
                default_value=_get_suffix,
                info="This suffix is appended (if not stated otherwise) to the SPRKKR "
                     "executable names."),
    V('dir', QString.I, is_optional=True, info='Directory, from which the executables will be runned. None mean use the default environment variable PATH mechanism')
  ], info="Configuration, that affects how the execubables are runned")

])

config = definition.create_object()
