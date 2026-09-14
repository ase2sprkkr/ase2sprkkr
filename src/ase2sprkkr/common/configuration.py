"""This module contains just a base class for both configuration values - :class:`Options<ase2sprkkr.common.options.Option>`
and configuration containers - :class:`Sections<ase2sprkkr.common.configuration_containers.Section>`.
"""

from contextlib import contextmanager
from typing import Generator, List, Literal, Optional, Tuple, Union
from .configuration_transaction import (
    ConfigurationTransaction,
    PostChangeHookError,
)
from .warnings import (
    DataValidityError,
    DataValidityErrors,
    DataValidityWarning,
    InvalidValuePolicy,
    ReportInvalidPolicy,
    RetainInvalidPolicy,
    ValidationReason,
    ValidationResult,
)
import warnings


UnknownMemberPolicy = Optional[Literal["find", "add", "ignore", "fail"]]


class Configuration:
    """The common base class for all configurations values and containers. I.e.
    for :class:`Options<ase2sprkkr.common.options.Option>` and :class:`Sections<ase2sprkkr.common.configuration_containers.Section>`.
    """

    def __init__(self, definition, container=None):
        """
        Create the object. Just sets the two properties from the parameters.

        Parameters
        ----------
        definition: ase2sprkkr.common.configuration_definitions.BaseDefinition
          Definition of this configuration object.

        container: ase2sprkkr.common.configuration_containers.ConfigurationContainer
          The container, that owns this configuration object.
        """

        self._definition = definition
        """
      The "definition" of the option or section. The definition determines
      the name(s), value type(s) etc... contained in the configuration object.
      Instance of :class:`ase2sprkkr.common.configuration_definitions.BaseDefinition`
      """

        self._container = container
        """
      The parent container. I.e. the container that holds this object (e.g.
      for a value it is the section that owns the value)
      Instance of :class:`ase2sprkkr.common.configuration_containers.ConfigurationContainer`
      """

    def _get_path(self, include_root=False):
        """Return the dot-delimited path to the item in the configuration tree.

        E.g. the ``ENERGY`` option in the ``CONTROL`` section has the path
        ``CONTROL.ENERGY``
        """
        name = self.name
        if self._container and (include_root or self._container._container):
            return f"{self._container._get_path()}.{name}"
        return name

    def _get_root_container(self):
        """Return the root object of the configuration tree.

        I.E. the object, that represents the whole configuration or problem-definition file
        """
        return self._container._get_root_container() if self._container else self

    def clear(self) -> None:
        """Transactionally clear this item and validate the resulting tree."""
        with self._mutation("set") as (transaction, _invalid):
            self.stage_clear(transaction)

    @property
    def name(self):
        """Return the name of the option/section. The name is defined by the definition
        of the object.

        Returns
        -------
        name: str
        The name of the object.
        """
        return self._definition.name

    def _as_dict(self, get):
        raise NotImplementedError()

    @staticmethod
    def as_dict_getter(only_changed: Union[bool, str] = "default", generated=False, copy=False):

        if only_changed == "default":

            def pick_only_changed(d):
                return not d.is_always_added
        elif only_changed == "basic":

            def pick_only_changed(d):
                return d.is_expert
        else:
            only_changed = bool(only_changed)

            def pick_only_changed(d):
                return only_changed

        def get(self):
            d = self._definition
            if d.is_generated and not generated:
                return None
            if only_changed == "explicit":
                v = self._unpack_value(self._value)
            elif pick_only_changed(d):
                v, c = self.value_and_changed()
                if not c:
                    return None
            else:
                v = self(all_values=True)
            if v is not None:
                if copy:
                    v = self._definition.copy_value(v, all_values=True)
            return v

        return get

    def as_dict(self, only_changed: Union[bool, str] = "basic", generated=False, copy=False, getter=None):
        """Return the value of self, in the case of container as a dictionary. To be redefined in the descendants.

        Parameters
        ----------
        only_changed
          Return only changed values, or all of them?
          If True, return only the values, that differ from the defaults.
          If False, return all the values.
          The default value ``basic`` means, return all non-expert values
          and all the changed expert values.
          ``explicit`` means just the values, that were explicitly set (even
          if they are the same as the default value)
        """
        if not getter:
            getter = self.as_dict_getter(only_changed, generated, copy)
        return self._as_dict(getter)

    to_dict = as_dict

    def show(self):
        """Print the configuration, as it will be saved into the configuration/problem definition file."""
        print(self.to_string())

    @property
    def info(self):
        return self._definition.info()

    @property
    def doc(self):
        try:
            return self._definition.description()
        except AttributeError as e:
            raise Exception("Cannot retrieve documentation") from e

    def help(self, verbose=False, show_hidden=False):
        if verbose is True:
            verbose = "all"
        elif verbose is False:
            verbose = True
        print(self._definition.description(verbose, show_hidden))

        global _help_warning_printed
        if not _help_warning_printed:
            import __main__ as main

            if verbose is True and not hasattr(main, "__file__"):  # I'm in repl
                print(
                    "\n You can use <Configuration>.help(True) for a more detailed description of the possible configuration options. Enjoy ASE2SPRKKR!\n"
                )
            _help_warning_printed = True

    def __repr__(self):
        d = self._definition
        out = d.configuration_type_name
        out = out + " " + d.name.upper()
        return out

    def check_for_errors(
        self, why: str = "save", print: bool = True
    ) -> List[warnings.WarningMessage]:
        """Return validation warnings for phase ``why`` and optionally print them."""
        with warnings.catch_warnings(record=True) as found:
            warnings.simplefilter("always", DataValidityWarning)
            self.validate(why, report_invalid="warn")
        if print:
            for warning in found:
                warnings.showwarning(
                    message=warning.message,
                    category=warning.category,
                    filename=warning.filename,
                    lineno=warning.lineno,
                )
        return found

    def validate(
        self,
        why: ValidationReason = "save",
        report_invalid: Optional[ReportInvalidPolicy] = None,
    ) -> List[ValidationResult]:
        """Validate this configuration.

        Parameters
        ----------
        why
          Validation phase: ``set``, ``parse`` or ``save``.

        report_invalid
          ``"raise"``, ``"warn"`` or ``"ignore"``. By default, errors
          warn after parsing and raise while setting or saving.

        """
        if why not in ("set", "parse", "save"):
            raise ValueError(f"Unknown validation reason: {why!r}")
        if report_invalid is None:
            report_invalid = "warn" if why == "parse" else "raise"
        results = ValidationResult.collect(self._validate, why)
        ValidationResult.report(results, report_invalid)
        return results

    @contextmanager
    def _mutation(
        self,
        why: ValidationReason,
        retain_invalid: Optional[RetainInvalidPolicy] = None,
        report_invalid: Optional[ReportInvalidPolicy] = None,
    ) -> Generator[
        Tuple[ConfigurationTransaction, InvalidValuePolicy], None, None
    ]:
        """Stage, validate and finish one policy-controlled mutation.

        ``why`` is the validation phase. ``retain_invalid`` and
        ``report_invalid`` have the values documented by
        :class:`InvalidValuePolicy`.
        """
        policy = InvalidValuePolicy(why, retain_invalid, report_invalid)
        hook_error = None

        try:
            with ConfigurationTransaction.use(self, policy) as transaction:
                policy = transaction.policy
                yield transaction, policy
                root = self._get_root_container()
                policy.add_semantic(
                    ValidationResult.collect(root._validate, why)
                )
                if policy.rollback:
                    transaction.abort()
        except PostChangeHookError as error:
            hook_error = error

        try:
            policy.report()
        except DataValidityError as validation_error:
            if hook_error is not None:
                raise DataValidityErrors(
                    "Configuration validation and post-change hooks failed",
                    (validation_error, hook_error),
                ) from None
            raise

        if hook_error is not None:
            raise hook_error

    def save_to_file(self, file, *, validate: Union[str, bool] = "save"):
        """Save the configuration to a file in a given format.

        This routine do some basic stuff and then call _save_to_file routine,
        that contains the implementation specific for the type of the
        configuration container/value.

        Parameters
        ----------
        file: str or file
          File to read the data from

        validate
          Validation mode. ``True`` and ``"save"`` perform strict save
          validation, ``False`` skips validation, and ``"warning"`` performs
          save validation without raising :class:`DataValidityError` results.
          ``"set"`` can be used to allow values that are only required when
          saving.
        """
        if validate is True:
            self.validate("save")
        elif validate is False:
            pass
        elif validate == "warning":
            self.validate("save", report_invalid="warn")
        else:
            self.validate(validate)

        if not hasattr(file, "write"):
            with open(file, "w") as file:
                out = self._save_to_file(file)
                file.flush()
        else:
            out = self._save_to_file(file)
            file.flush()
        return out

    def to_string(self, *, validate: Union[str, bool] = "warning"):
        """
        Return the configuration (problem definition) in a string.

        Parameters
        ----------
        Validate. How to validate before retrieving. See the method
        validate.
        Default 'warning' means the same as 'save', but only throw a warning.

        Returns
        -------
        configuration:str
        The configuration, as it should be saved in a configuration/problem definition file.
        """
        from io import StringIO

        s = StringIO()
        self.save_to_file(s, validate=validate)
        return s.getvalue()

_help_warning_printed = False
