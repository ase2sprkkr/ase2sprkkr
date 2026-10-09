import asyncio
import os
import sys
import tempfile
import threading
import time

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ..process_output_reader import ProcessOutputParser  # NOQA: E402


class _ProcessParser(ProcessOutputParser):
    def __init__(self):
        super().__init__()
        self.stop_thread = None

    async def read_output(self, stdout, *args):
        await stdout.read()

    async def _stop_process(self):
        self.stop_thread = threading.get_ident()
        await super()._stop_process()

    def result(self, output, error, wait):
        return wait


class TestProcessOutputParser(TestCase):
    def test_stop_the_process_runs_on_owning_event_loop(self):
        parser = _ProcessParser()
        errors = []

        def run():
            try:
                parser._run(
                    [sys.executable, "-c", "import time; time.sleep(30)"],
                    None,
                )
            except BaseException as exc:
                errors.append(exc)

        thread = threading.Thread(target=run, daemon=True)
        thread.start()

        deadline = time.monotonic() + 5
        while parser.proc is None and time.monotonic() < deadline:
            time.sleep(0.01)
        self.assertTrue(parser.proc is not None)

        caller_thread = threading.get_ident()
        parser.stop_the_process()
        thread.join(5)

        self.assertFalse(thread.is_alive())
        self.assertEqual(errors, [])
        self.assertEqual(parser.stop_thread, thread.ident)
        self.assertNotEqual(parser.stop_thread, caller_thread)

    def test_stop_before_run_prevents_process_start_and_restores_cwd(self):
        parser = _ProcessParser()
        original_cwd = os.getcwd()
        parser.stop_the_process()

        with tempfile.TemporaryDirectory() as directory:
            self.assertRaises(
                asyncio.CancelledError,
                lambda: parser._run(
                    [sys.executable, "-c", "pass"],
                    None,
                    directory=directory,
                ),
            )

        self.assertIsNone(parser.proc)
        self.assertEqual(os.getcwd(), original_cwd)
