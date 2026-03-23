""" The reader for JXC task, providing exchange couplings, exchange tensors, and Dzyaloshinski-Moriya vectors. """
import warnings
import re

from ..task_result import TaskResult, KkrProcessRunner, OutputFileResultValue
from .default import SprKkrOutputReader
from ...common.decorators import cached_property
from ...bindings.uppasd import write_pos_file, write_mom_file
from ...common.process_output_reader import readline_until


class JxcResult(TaskResult):
    """ JXC task provides access to the generated JXC output file.
    It provides access to exchange couplings, exchange tensors, and Dzyaloshinski-Moriya vectors,
    using the :py:attr:`~jxc`, :py:attr:`~dij` and :py:attr:`~dmi` properties.
    Not all of them may be present in the output, depending on the settings of the JXC task.
    """

    @cached_property
    def dmi_filename(self):
        """ Dzyaloshinski-Moriya vectors file name """
        try:
            return self.path_to('dmi')
        except KeyError:
            raise AttributeError('Dzyaloshinski-Moriya vectors file not found in output')

    @cached_property
    def dmi(self):
        """ Dzyaloshinski-Moriya vectors. """
        return self.output_values['dmi']()

    @cached_property
    def dij_filename(self):
        """ Exchange tensors file name """
        try:
            return self.path_to('dij')
        except KeyError:
            raise AttributeError('Exchange tensors file not found in output')

    @cached_property
    def dij(self):
        """ Exchange tensors. """
        return self.output_values['dij']()

    @cached_property
    def jxc_filename(self):
        """ Exchange couplings file name """
        try:
            return self.path_to('jxc')
        except KeyError:
            raise AttributeError('Exchange couplings file not found in output')

    @cached_property
    def jxc(self):
        """ Exchange couplings. """
        return self.output_values['jxc']()

    @property
    def output_values(self):
        files = {
            'jxc': 'Exchange couplings',
            'dij': 'Exchange tensors',
            'dmi': 'Dzyaloshinski-Moriya vectors',
        }

        return {
            key: OutputFileResultValue(name, key, getattr(self, f'{key}_filename'))
            for key, name in files.items()
            if key in self.files
        }

    def write_uppasd_files(self, directory=None):
        """ Write the computed exchange couplings, exchange tensors and Dzyaloshinski-Moriya vectors
        to UppASD input files in the specified directory. """
        output_file = None
        for key in ['jxc', 'dmi']:
            if key in self.files:
                output_file = self.output_values[key]()
                output_file.write_uppasd_file(directory=directory)

        write_pos_file(self.atoms, jxc_ouput_file=output_file, directory=directory)
        write_mom_file(self.atoms, jxc_ouput_file=output_file, directory=directory)


class JxcOutputReader(SprKkrOutputReader):

     async def read_output(self, stdout, result):
        await self.read_commons(stdout, result)

        pattern = re.compile(r'\s*results written to file:(.+_(XCPLTEN_Jij|XCPLTEN_Dij|DMIVEC_Dij)\.dat)\s*')
        file_types = {'XCPLTEN_Jij': 'jxc', 'XCPLTEN_Dij': 'dij', 'DMIVEC_Dij': 'dmi'}

        while True:
            line = await readline_until(stdout, lambda line: line.startswith(b'          results written to file:'))
            if line is None:
                break
            match = pattern.match(line)
            if match:
                result.files[file_types[match.group(2)]] = match.group(1)
            else:
                breakpoint()
                warnings.warn(f'Unexpected line in JXC output: {line}')

class JxcProcessRunner(KkrProcessRunner):

    result_class = JxcResult
    reader_class = JxcOutputReader
