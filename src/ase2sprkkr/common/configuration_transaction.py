"""Transactional staging and committing of configuration changes."""

from __future__ import annotations

from types import TracebackType
from typing import Any, Callable, Dict, Optional, Sequence, Type, TYPE_CHECKING

from .backward_compatibility import ExceptionGroup
from .warnings import InvalidValuePolicy

if TYPE_CHECKING:
    from .configuration import Configuration


class PostChangeHookError(RuntimeError, ExceptionGroup):
    """One or more post-change hooks failed after a successful commit."""

    def __init__(self, message: str, exceptions: Sequence[Exception]) -> None:
        """Initialize the exception-group base explicitly for Python 3.8."""
        ExceptionGroup.__init__(self, message, exceptions)

    def derive(self, exceptions: Sequence[Exception]) -> ExceptionGroup:
        """Create the same group type for the selected ``exceptions``."""
        return type(self)(self.message, exceptions)


class ConfigurationTransaction:
    """Context manager for temporary configuration mutations.

    Staged objects hold their proposed state directly, so ordinary container
    traversal sees the complete proposal. A normal context exit commits all
    staged objects; an exceptional exit rolls them back in reverse order.

    Mutable values are the real objects. Their in-place mutations deliberately
    cannot be rolled back by this transaction.
    """

    class _Boundary:
        """Change-stack marker and explicit savepoint context."""

        def __init__(self, transaction: "ConfigurationTransaction") -> None:
            """Create an inactive boundary owned by ``transaction``."""
            self._transaction = transaction
            self._active = False
            self._rollback_requested = False

        def __enter__(self) -> "ConfigurationTransaction._Boundary":
            if self._active:
                raise RuntimeError("The savepoint is already active")
            self._rollback_requested = False
            self._active = True
            transaction = self._transaction
            if not transaction._changes:
                transaction._container._active_transaction = transaction
            transaction._changes.append(self)
            return self

        def __exit__(
            self,
            exception_type: Optional[Type[BaseException]],
            exception: Optional[BaseException],
            traceback: Optional[TracebackType],
        ) -> bool:
            return self._transaction._leave(
                self,
                exception_type is not None or self._rollback_requested,
            )

        def rollback(self):
            """Rollback changes registered since the savepoint was created."""
            if not self._active:
                raise RuntimeError("The savepoint is not active")
            self._transaction._rollback_after(self)
            self._rollback_requested = True

    def __init__(
        self,
        configuration: "Configuration",
        policy: Optional[InvalidValuePolicy] = None,
    ) -> None:
        """Create a transaction for ``configuration`` and its root tree.

        ``policy`` is shared by every staged value in a high-level mutation.
        When omitted, the transaction uses a strict policy for direct,
        low-level staging.
        """
        self._container = configuration._get_root_container()
        if getattr(self._container, "_running_post_change_hooks", False):
            raise RuntimeError(
                "A configuration cannot be changed from its post-change hooks"
            )
        self.policy = policy or InvalidValuePolicy(
            "set", retain_invalid="none", report_invalid="raise"
        )
        self._defer_validation = policy is not None
        self._changes = []

    @classmethod
    def current(
        cls, configuration: "Configuration"
    ) -> Optional["ConfigurationTransaction"]:
        """Return the active transaction for ``configuration``, if any."""
        root = configuration._get_root_container()
        return getattr(root, "_active_transaction", None)

    @classmethod
    def use(
        cls,
        configuration: "Configuration",
        policy: Optional[InvalidValuePolicy] = None,
    ) -> "ConfigurationTransaction":
        """Return a context for the transaction applicable to ``configuration``.

        An active transaction is entered as a nested savepoint; otherwise a
        new transaction is returned. A failed nested operation rolls back only
        its own changes unless its exception escapes the outer transaction.
        ``policy`` is used only when a new transaction is created; nested
        operations share the policy of the active transaction.
        """
        transaction = cls.current(configuration)
        if transaction is not None:
            return transaction
        return cls(configuration, policy=policy)

    def _finish_stage(self) -> None:
        """Report a standalone stage immediately.

        Transactions created directly are strict low-level transactions. A
        policy supplied by :meth:`Configuration._mutation` is reported only
        after the complete mutation, so its results can be aggregated.
        """
        if self._defer_validation:
            return
        try:
            self.policy.report()
        finally:
            self.policy.results.clear()

    def __enter__(self) -> "ConfigurationTransaction":
        self._Boundary(self).__enter__()
        return self

    def __exit__(
        self,
        exception_type: Optional[Type[BaseException]],
        exception: Optional[BaseException],
        traceback: Optional[TracebackType],
    ) -> bool:
        boundary = self._current_boundary()
        return self._leave(
            boundary,
            exception_type is not None or boundary._rollback_requested,
        )

    def _leave(
        self,
        boundary: "ConfigurationTransaction._Boundary",
        rollback: bool,
    ) -> bool:
        """Close ``boundary``, rolling it back when ``rollback`` is true."""
        self._require_active()
        outermost = self._changes[0] is boundary
        boundary._active = False
        hooks = None
        try:
            if rollback:
                self._rollback_after(boundary)
            if outermost:
                if not rollback:
                    hooks = self._commit()
        finally:
            if outermost:
                self._changes.clear()
                self._container.__dict__.pop("_active_transaction", None)
        if hooks:
            self._run_hooks(hooks)
        return False

    def _require_active(self) -> None:
        if not self._changes:
            raise RuntimeError(
                "A configuration transaction can only be used inside its "
                "active context manager"
            )

    def push(self, change: Any) -> Any:
        """Push a staged change and return it; it must implement rollback."""
        try:
            self._require_active()
        except Exception:
            change.rollback()
            raise
        self._changes.append(change)
        return change

    def savepoint(self) -> "ConfigurationTransaction._Boundary":
        """Rollback newly registered changes if the enclosed operation fails."""
        self._require_active()
        return self._Boundary(self)

    def abort(self) -> None:
        """Roll back changes in the current transaction boundary."""
        self._current_boundary().rollback()

    def _current_boundary(self) -> "ConfigurationTransaction._Boundary":
        for change in reversed(self._changes):
            if isinstance(change, self._Boundary) and change._active:
                return change
        raise RuntimeError("The transaction has no active boundary")

    def _rollback_after(
        self, boundary: "ConfigurationTransaction._Boundary"
    ) -> None:
        """Rollback changes above ``boundary``, leaving it active."""
        while self._changes and self._changes[-1] is not boundary:
            change = self._changes.pop()
            if isinstance(change, self._Boundary):
                if change._active:
                    raise RuntimeError(
                        "A nested transaction boundary is still active"
                    )
            else:
                change.rollback()
        if not self._changes:
            raise RuntimeError("The transaction boundary is no longer active")

    def _commit(self) -> Dict[Any, Callable[[Any], None]]:
        """Commit staged objects and collect their post-change hooks."""
        hooks = {}
        while self._changes:
            changes = self._changes
            self._changes = []
            for change in changes:
                if isinstance(change, self._Boundary):
                    continue
                hook = change.commit()
                if hook:
                    option, callback = hook
                    hooks[option] = callback
        return hooks

    def _run_hooks(self, hooks: Dict[Any, Callable[[Any], None]]) -> None:
        """Run registered ``hooks`` after closing the transaction."""
        errors = []
        self._container._running_post_change_hooks = True
        try:
            for option, hook in hooks.items():
                try:
                    hook(option)
                except Exception as error:
                    errors.append(error)
        finally:
            self._container.__dict__.pop("_running_post_change_hooks", None)
        if errors:
            raise PostChangeHookError(
                "Configuration was already committed, but post-change hooks failed",
                errors,
            )

    def __repr__(self) -> str:
        return f"<Configuration transaction for {self._container}>"
