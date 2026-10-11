import tempfile

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ..process_output_reader import ProcessOutputParser  # NOQA: E402


class _OutputParser(ProcessOutputParser):
    async def read_output(self, stdout, *args):
        out = b""
        while line := await stdout.readline():
            out += line
        return out

    async def read_error(self, stderr, *args):
        out = b""
        while line := await stderr.readline():
            out += line
        return out


class TestProcessOutputParser(TestCase):
    def test_read_output_file(self):
        parser = _OutputParser()
        with tempfile.NamedTemporaryFile() as output, tempfile.NamedTemporaryFile() as error:
            output.write(b"one\ntwo\n")
            output.flush()
            error.write(b"warning\n")
            error.flush()

            out, err, return_code = parser.read_from_file(
                output.name, error.name
            )

        self.assertEqual(out, b"one\ntwo\n")
        self.assertEqual(err, b"warning\n")
        self.assertEqual(return_code, 0)
