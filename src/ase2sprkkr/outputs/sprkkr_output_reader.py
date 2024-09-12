from ..common.process_output_reader import ProcessOutputReader, readline_until
import os
import datetime
import re


class SprKkrOutputReader(ProcessOutputReader):

    async def read_commons(self, stdout, result):
        out = await self.parse_files(stdout, result)
        return out

    async def parse_files(self, stdout, result):
        version = await readline_until(stdout, lambda line: b'VERSION' in line, can_end=False)
        version = version.split()
        result.program_info = {'version' : version[3],
                               'executable' : version[1] }
        started = await readline_until(stdout, lambda line: b'programm execution' in line, can_end=False)
        started = re.sub('[a-z]','', started).strip()
        result.program_info['start_time']=datetime.datetime.strptime(started, "%d/%m/%Y %H:%M:%S")

        await readline_until(stdout, lambda line: line == b' fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff\n', can_end=False)
        await stdout.readline()
        out=await stdout.readline()
        if out.decode('utf8').strip() != 'FILES:':
            for i in range(10):
              print(await stdout.readline())
            raise ValueError(f"Unexpected line: {out}")
        await stdout.readline()
        line = (await stdout.readline())
        line=line.decode('utf8').strip()
        while line:
            file=line.split(':')
            result.files[file[0].strip()]=file[1].split(')',1)[1].strip()
            line = await stdout.readline()
            line=line.decode('utf8').strip()

        if not 'input' in result.files:
            if hasattr(stdout, 'file') and hasattr(stdout.file, 'name'):
                filename = stdout.file.name
                if filename.endswith('.out') and os.path.exists(filename[:-4] + '.inp'):
                    result.files['input'] = filename[:-4] + '.inp'
