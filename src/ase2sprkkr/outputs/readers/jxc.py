""" The reader for JXC task, providing exchange couplings, exchange tensors, and Dzyaloshinski-Moriya vectors. """
import warnings
import re

from ..task_result import TaskResult, KkrOutputReader, OutputFileResultValue
from ..sprkkr_output_reader import SprKkrOutputParser
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
            key: OutputFileResultValue(name, 'jxc', getattr(self, f'{key}_filename'))
            for key, name in files.items()
            if key in self.files
        }

    def write_uppasd_files(self, directory=None, jxc_file_name='jfile.dat', dmi_file_name='dmfile.dat', pos_file_name='posfile.dat', mom_file_name='momfile.dat'):
        """ Write the computed exchange couplings, exchange tensors and Dzyaloshinski-Moriya vectors
        to UppASD input files in the specified directory. """
        output_file = None
        for key  in ['jxc', 'dmi']:
            if key in self.files:
                output_file = self.output_values[key]()
                file_name = locals()[key + '_file_name']
                output_file.write_uppasd_file(file_name=file_name, directory=directory)

        write_pos_file(self.atoms, pos_file_name, jxc_ouput_file=output_file, directory=directory)
        write_mom_file(self.atoms, mom_file_name, jxc_ouput_file=output_file, directory=directory)


class JxcOutputParser(SprKkrOutputParser):

     async def read_output(self, stdout, result):
        await self.read_commons(stdout, result)

        pattern = re.compile(rb'\s*(?:(?:' + rb')|(?:'.join([
            rb'(?P<file>results written to file:(?P<file_name>.+_(?P<file_type>XCPLTEN_Jij|XCPLTEN_Dij|DMIVEC_Dij)\.dat|.*))',
            rb'Curie temperature within mean field approximation  T_C =\s*(?P<curie_temp>[0-9.]+) K',
        ])+rb'))\s*')
        file_types = {b'XCPLTEN_Jij': 'jxc', b'XCPLTEN_Dij': 'dij', b'DMIVEC_Dij': 'dmi'}

        match = None
        def important(line):
            nonlocal match
            match = pattern.match(line)
            return match

        result.mean_field_curie_temperature = None
        while True:
            line = await readline_until(stdout, important)
            if line is None:
                break
            if match.group("file"):
                file_type = file_types.get(match.group("file_type"), None)
                if file_type is not None:
                    result.files[file_type] = match.group("file_name").decode('utf8')
                else:
                    warnings.warn(f'Unexpected file type in JXC output: {line}')
            if match.group("curie_temp"):
                result.mean_field_curie_temperature = float(match.group("curie_temp"))

        if result.mean_field_curie_temperature is None:
            warnings.warn('Curie temperature not found in JXC output')


class JxcOutputReader(KkrOutputReader):

    result_class = JxcResult
    parser_class = JxcOutputParser
