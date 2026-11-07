# Extraction validation

Source baselines: ASE2SPRKKR `b711659b794e911a2356fca8d53d7c5a27775891`,
Guy4ASE `b1d18e1fa64ced0314e434f871d66869cc5efea7`.

| Check | Result |
| --- | --- |
| KKRunner complete default suite | 678 passed, 3 optional Docker cases skipped |
| Real SSH reconnect (immediately after submit; after first stdout) and PSI/J Slurm | 3 passed |
| ASE2SPRKKR execution and sprkkr groups | 362 passed, 4 skipped |
| Additional preparation/resource tests including transient-path rejection | 20 passed |
| ASE2SPRKKR input-parameter suite | 1126 passed |
| ASE2SPRKKR output readers | 76 passed |
| ASE2SPRKKR end-user remote CLI | 8 passed |
| Guy4ASE full suite after integration | 789 passed, 24 failed, 19 skipped |
| Guy4ASE full baseline suite | 787 passed, 26 failed, 19 skipped |

Every remaining GUI failure ID was also present in the baseline. Examples are
missing `choose_directory` arguments in existing widget construction, an existing
SPLEED/ARPES default-schema assertion, and existing recompute keyword expectations.
No new GUI failure IDs were introduced. The two former recovery-boundary failures
now pass. These unrelated GUI defects were not changed to make the suite green.

The broad ASE2SPRKKR invocation and original execution baseline exceeded their
initial time limits. Relevant current groups were subsequently run separately;
the table reports completed groups, not a claim that the full project suite is
clean. Native examples and optional external-service cases are separate checks.

A pure Python wheel was built and installed into a fresh virtual environment
containing only KKRunner and platformdirs. Both `ase` and `ase2sprkkr` were absent.
Raw local execution, persistent execution with a callback supplied after restart,
and the remote CLI/service via a separate interpreter all passed there.
The dependency wheel was reconstructed from the locally installed platformdirs
distribution to keep this installation test offline.

Real SPR-KKR 9.7 executables were used for a small Al SCF run (2 iterations,
NKTAB=5, NE=5) and a subsequent DOS run (NKTAB=5, NE=5). Both returned 0,
finished their persistent runs and were loaded again with TaskResult.from_run.
Both recorded task identity, schema identity and verified output artifacts.

The SSH/Slurm image installs no ASE/ASE2SPRKKR. The SSH case explicitly checks
that neither application package is present on the target. The image receives
source through Docker build context rather than requiring host-visible binds.

Python 3.14 was used for host tests; the Docker integration uses Python 3.12.
All production modules also parse with Python 3.9's syntax grammar. This is not
a claim of a complete Python 3.9 runtime test matrix.

Ruff passes for all runner production modules and the changed execution/result
adapters. The pre-existing unused `canonical` assignment in the application
input-definition registry is outside that clean subset. git diff --check passes.

## Reproduction

```sh
# from the kkrunner repository
python -m pip install -e '.[test]'
python -m pytest -q
KKRUNNER_RUN_DOCKER_TESTS=1 python -m pytest -q tests/test_ssh_standalone.py tests/test_psij_slurm.py

# from ASE2SPRKKR after installing the sibling runner
DO_NOT_RUN_SPRKKR=1 python -m pytest -q src/ase2sprkkr/input_parameters/test
DO_NOT_RUN_SPRKKR=1 python -m pytest -q src/ase2sprkkr/outputs/test
DO_NOT_RUN_SPRKKR=1 python -m pytest -q src/ase2sprkkr/execution/test src/ase2sprkkr/sprkkr/test \
  --ignore=src/ase2sprkkr/execution/test/test_remote_ssh.py \
  --ignore=src/ase2sprkkr/execution/test/test_remote_psij_slurm.py

# from Guy4ASE
QT_QPA_PLATFORM=offscreen DO_NOT_RUN_SPRKKR=1 python -m pytest -q tests
```

Execution integration tests require subprocesses and Unix sockets. An environment
that prohibits those operations cannot exercise supervisor/reconnection behavior.
