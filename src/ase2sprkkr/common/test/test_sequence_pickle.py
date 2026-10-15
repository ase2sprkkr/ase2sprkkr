import multiprocessing
import pickle

from ase2sprkkr.common.grammar_types import Sequence


def _check_named_sequence_in_spawned_process(sequence, value):
    assert value.VERSION == 9
    assert value.DATE == "2019-01-18"
    assert isinstance(value, sequence.value_type)

    rebuilt = sequence.value_constructor(10, "2026-10-10")
    assert rebuilt.VERSION == 10
    assert rebuilt.DATE == "2026-10-10"


def test_named_sequence_and_value_are_pickleable():
    sequence = Sequence(int, str, names=("VERSION", "DATE"))
    value = sequence.value_constructor(9, "2019-01-18")

    restored_sequence = pickle.loads(pickle.dumps(sequence))
    restored_value = pickle.loads(pickle.dumps(value))

    assert restored_value == value
    assert isinstance(restored_value, restored_sequence.value_type)


def test_named_sequence_and_value_can_be_sent_to_spawned_process():
    sequence = Sequence(int, str, names=("VERSION", "DATE"))
    value = sequence.value_constructor(9, "2019-01-18")

    context = multiprocessing.get_context("spawn")
    process = context.Process(
        target=_check_named_sequence_in_spawned_process,
        args=(sequence, value),
    )
    process.start()
    process.join(10)
    if process.is_alive():
        process.terminate()
        process.join()
        raise AssertionError("Spawned process did not finish")
    assert process.exitcode == 0
