from ..process_output_reader import ProcessOutputParser


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


def test_read_output_file(tmp_path):
    output = tmp_path / "stdout"
    error = tmp_path / "stderr"
    output.write_bytes(b"one\ntwo\n")
    error.write_bytes(b"warning\n")

    out, err, return_code = _OutputParser().read_from_file(output, error)

    assert out == b"one\ntwo\n"
    assert err == b"warning\n"
    assert return_code == 0
