from ..common.process_output_reader import ProcessOutputParser, readline_until
from ..common.file_utils import filename_from_file
import os
import datetime
import re


class SprKkrOutputParser(ProcessOutputParser):
    _common_output_line = re.compile(rb"VERSION|programm execution|^\s*from file:|^\s*FILES:\s*$")
    _sfn_file = re.compile(r"^\s*from file:\s*(.*?)\s*$")

    async def read_commons(self, stdout, result):
        out = await self.parse_files(stdout, result)
        return out

    async def parse_files(self, stdout, result):
        # Read version info
        try:
            result.program_info = {}
            while True:
                line = await readline_until(stdout, self._common_output_line.search, can_end=False)

                if "version" not in result.program_info and "VERSION" in line:
                    version = line.split()
                    result.program_info.update(version=version[3], executable=version[1])
                elif "programm execution" in line:
                    started = re.sub("[a-z]", "", line).strip()
                    result.program_info["start_time"] = datetime.datetime.strptime(started, "%d/%m/%Y %H:%M:%S")
                elif match := self._sfn_file.match(line):
                    result.files.add_file("SFN", match.group(1), "sfn")
                elif line.strip() == "FILES:":
                    break

            await stdout.readline()  # skip empty line

            # Parse files
            line = (await stdout.readline()).decode("utf8").strip()
            while line:
                if ":" in line:
                    name, rhs = line.split(":", 1)
                    name = name.strip()
                    rhs = rhs.strip()

                    # Remove (number) if present
                    if rhs.startswith("(") and ")" in rhs:
                        rhs = rhs.split(")", 1)[1].strip()

                    result.files.add_file(name, rhs)

                line = (await stdout.readline()).decode("utf8").strip()

            # Fallback for input file
            if "input" not in result.files:
                filename = filename_from_file(stdout, None)
                if filename:
                    if filename.endswith(".out") and os.path.exists(filename[:-4] + ".inp"):
                        result.files.add_file("input", filename[:-4] + ".inp")

        except EOFError as e:
            raise EOFError("Unexpected end of output -- the program exited prematurely.") from e
