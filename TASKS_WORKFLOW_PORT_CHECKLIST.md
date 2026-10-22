# Port of `tasks` onto `workflow`

> **Modernized on 2026-10-22:** This branch has been rebuilt from the
> current `workflow` tip (`dfd21d18bac1`).
> A three-way comparison with the previous `tasks-workflow` tip
> (`b15e2cef84dd`) found 208 files with identical task-related changes,
> 53 with independent newer changes on `workflow`, and one task-only
> regression file (`test_persistent_option.py`). The newer `workflow`
> implementations are retained for all 53 divergent files to avoid rolling
> back runner, registry, validation, or SPC/NeXus fixes.
> This updates the **base**, not the completion status of the original
> `tasks` parity audit; outstanding checks below remain outstanding.
> Validate the new test and full CI before integrating.

This file tracks **all** functional changes to migrate. Completion requires
comparison against the original `tasks` branch, not only a passing CI job.

## Architecture and common parsing

- [x] Separate original `tasks` from `tasks-workflow` (modern base).
- [x] Transfer preliminary independent task definitions with shared schemas.
- [ ] Verify all migrated task factories construct, roundtrip, validate and pickle.
- [ ] Support **nonempty bracketed string sets** in common input `Mixed` via
      `SetOf(require_brackets=True)`; verify CI. Remove parser-only `VariantGenerator`.
- [ ] Check other native `Mixed` array shapes (numbers, booleans, quoted strings,
      braces and parentheses), absence of ambiguous unbracketed parsing, and empty sets.
- [ ] Keep exact IPD/task namespace rules and schema discovery lazy where intended.

## Per-task files: *declarative* task ownership

- [ ] Support file specifications with parameter-dependent conditions; preserve
      input/output/cache distinctions, wildcard templates, directory resolution.
- [ ] Implement KKRSPEC auxiliary input `input.inp` if `INPVER==0` and
      `struc.inp` if `STRVER==0`.
- [ ] Implement KKRSPEC `{dataset}_SPEC.out` and
      `{dataset}_data.spc` versus `_data_pol1.spc`/`_data_pol2.spc`,
      conditional on native task and polarization.
- [ ] Apply this consistently to ARPES/AIPES/BAND/SPLEED/CLPES.
- [ ] Audit **every** task's required inputs, declared outputs, optional cache,
      filename patterns and parameter-dependent conditions against `tasks`.
- [ ] Check overlap handling, empty outputs, and wildcard discovery behavior
      against the current task/files infrastructure.

## Existing definitions: no silent omissions

- [ ] SCF: transfer and test additional range, mixing, grid, energy and
      FULLPOT/TCHEBY validations.
- [ ] DOS: transfer and test ImE, grids, NE and finite energy bounds.
- [ ] BSF: diff validations, modes, k-paths and cache/output manifest.
- [ ] ARPES: diff full SPEC_INPUT validation and file behavior.
- [ ] PHAGEN/XAS/JXC/DMI and other *already present* definitions: line-by-line
      compare task-specific validators, defaults and declared files.
- [ ] Ensure native mode flags and physical validation are independent of
      requested task identity, except task-specific defaults where intended.

## Remaining native selector families

- [ ] PLOTPRS selectors PLOTPOT, PLOTRHO, PLOTSFN: shared IPD,
      distinct tasks, outputs and FULLPOT restrictions.
- [ ] PLOT2DPR selectors PLOT2DPOT, PLOT2DRHO: shared IPD,
      distinct tasks, geometry and FULLPOT behavior.
- [ ] AIPES/BAND/SPLEED: compare modern task defaults, validation and files;
      do not duplicate ARPES schema.
- [ ] Move appropriate `tasks` test cases into the modern test suite
      without bringing back deprecated registries.

## Tests, audit and completion

- [ ] Adapt old tests for the modern API (`p.task_definition`, not IPD
      execution metadata or `p.is_mpi()`).
- [ ] Test that each registered TaskDefinition points to an appropriate IPD,
      and task settings roundtrip and resolve to the correct task.
- [ ] Transfer/replace missing task-specific regression tests for
      AIPES, BAND, SPLEED, SCF, DOS, XAS, PHAGEN, plot families,
      file manifests and audit tests.
- [ ] Fix remaining substantive test failures: APS invalid shell,
      WFPLOT MJ warnings, CHIXAS physically required input, etc.
- [ ] Run full Python CI matrix, CLI execution-integration and Guy4ASE
      compatibility checks, inspect actual statuses.
- [ ] Independently audit port coverage vs the original branch,
      review exceptions and regressions, and clean commit history.
- [ ] Only merge after explicit user request and thorough audit.

**Status:** in progress. `tasks` and `workflow` remain unchanged.
