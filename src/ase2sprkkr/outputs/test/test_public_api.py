def test_persistent_result_helpers_are_public():
    from ase2sprkkr.outputs import (
        PersistentResultMismatchError,
        PersistentResultUnavailableError,
        TaskResult,
        UnrecoverablePersistentResultError,
        attach_persistent_run,
        check_persistent_result,
        load_persistent_result,
    )
    from ase2sprkkr.outputs.persistent import (
        PersistentResultMismatchError as implementation_mismatch_error,
        PersistentResultUnavailableError as implementation_error,
        UnrecoverablePersistentResultError as implementation_unrecoverable_error,
        attach_persistent_run as implementation_attach,
        check_persistent_result as implementation_check,
        load_persistent_result as implementation_load,
    )
    from ase2sprkkr.outputs.task_result import TaskResult as implementation_result

    assert PersistentResultMismatchError is implementation_mismatch_error
    assert PersistentResultUnavailableError is implementation_error
    assert UnrecoverablePersistentResultError is implementation_unrecoverable_error
    assert issubclass(PersistentResultMismatchError, PersistentResultUnavailableError)
    assert issubclass(UnrecoverablePersistentResultError, PersistentResultMismatchError)
    assert TaskResult is implementation_result
    assert attach_persistent_run is implementation_attach
    assert check_persistent_result is implementation_check
    assert load_persistent_result is implementation_load
