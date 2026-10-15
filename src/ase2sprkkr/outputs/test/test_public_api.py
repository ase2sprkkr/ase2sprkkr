def test_persistent_result_helpers_are_public():
    from ase2sprkkr.outputs import (
        PersistentResultUnavailableError,
        TaskResult,
        attach_persistent_run,
        load_persistent_result,
    )
    from ase2sprkkr.outputs.persistent import (
        PersistentResultUnavailableError as implementation_error,
        attach_persistent_run as implementation_attach,
        load_persistent_result as implementation_load,
    )
    from ase2sprkkr.outputs.task_result import TaskResult as implementation_result

    assert PersistentResultUnavailableError is implementation_error
    assert TaskResult is implementation_result
    assert attach_persistent_run is implementation_attach
    assert load_persistent_result is implementation_load
