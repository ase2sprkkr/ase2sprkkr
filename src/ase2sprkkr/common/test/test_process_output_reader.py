import io

from ..process_output_reader import AsyncioFileReader, ProcessOutputParser


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


def test_read_output_file_does_not_close_caller_streams():
    output = io.BytesIO(b"output\n")
    error = io.BytesIO(b"warning\n")

    parsed = _OutputParser().read_from_file(output, error)

    assert parsed == (b"output\n", b"warning\n", 0)
    assert not output.closed and not error.closed
    output.seek(0)
    error.seek(0)
    assert output.read() == b"output\n"
    assert error.read() == b"warning\n"


def test_reader_only_closes_streams_that_it_opens(tmp_path):
    source = tmp_path / "output"
    source.write_bytes(b"data")

    with AsyncioFileReader(source) as reader:
        owned = reader.file
        assert not owned.closed
    assert owned.closed

    with source.open("rb") as borrowed:
        with AsyncioFileReader(borrowed) as reader:
            assert reader.file is borrowed
            assert borrowed.read() == b"data"
        assert reader.file is None
        assert not borrowed.closed
