# Schválený výběr pro oficiální dokumentaci

**Stav: schváleno uživatelem 8. 11. 2026 s úpravou instalačního příkladu.**

Pevná relativní cesta ke checkoutu byla nahrazena citovaným placeholderem
`"<path-to-kkrunner>"`; příklad nevyžaduje sousední repozitáře. Nic nebylo publikováno.

Historické AI review, designové snapshoty a validační protokoly jsou oddělené
od oficiální nápovědy. Jejich původní nálezy, počty testů a verdikty nejsou
považovány za současný stav. Níže je výběr, který navrhuji ponechat jako
udržované návody. Citace jsou přesné výňatky z připravených stránek; odkazy
vedou na celé znění v pracovním stromu.

## 1. Rozdělení odpovědností

Hranice knihoven z původní architektury; současné třídy a žádná kompatibilní jména.

[Celá stránka](/home/logik/skola/kkr/kkrunner/docs/architecture.md)

> KKRunner owns ExecutionSpec/ExecutionOutcome, TaskFile, ExecutorEndpoint,
> RunState/FailurePhase, Execution/WorkingDirectory, RunRecord/RunStore,
> FileEntry/FileMapping/FileManifest, FileState/FileVersion/FileFingerprint,
> ExecutionResources and ResourceOwnedRunner. Public exports live in `kkrunner`.
> `RunManager` tracks runs in an explicit RunStore without requiring a project.
> `ProjectExecutionManager` and ProjectRepository organize runs in projects.
> Both managers use the same lifecycle. Each public class has one name that
> describes its role.

## 2. Úplný manifest a rezervace

R03, R06, R11: pravidlo deklarací místo starého seznamu chyb či rezervace celého adresáře.

[Celá stránka](/home/logik/skola/kkr/kkrunner/docs/files.rst)

> The manifest is the complete declaration of what a run may read or write.
> KKRunner reserves these paths before staging or opening outputs. It does not
> reserve a whole directory merely because two calculations share it.

## 3. Povinné a volitelné výstupy

R12: rozlišit chybějící artefakt a chybu přenosu, aby uživatel zvolil správný postup.

[Celá stránka](/home/logik/skola/kkr/kkrunner/docs/files.rst)

> Successful durable completion captures file versions and SHA-256 digests.
> Historical applications should verify captured files before interpreting them.
> :class:`~kkrunner.HistoricalFileMismatchError` indicates that historical content
> was replaced. A missing required local output is an ``ARTIFACTS`` failure;
> a failed remote transfer is a ``DOWNLOAD`` failure. Optional missing outputs
> are recorded as skipped.

## 4. Obnova a oddělené výsledky

R01 a opravy výsledkového typu: raw outcome není vědecký výsledek, chybějící procesor není úspěšné zpracování.

[Celá stránka](/home/logik/skola/kkr/kkrunner/docs/persistent-runs.rst)

> ``process_result`` consumes an :class:`~kkrunner.ExecutionOutcome` and returns
> an application result. It is separate from the completion-notification callback.
> ``processing_required=True`` allows supplying the processor after reconnect.
> If the processor is unavailable or another client owns/performed processing,
> ``wait()`` raises :class:`~kkrunner.ProcessingPending`; it must not substitute a
> raw outcome for the application's result type.

## 5. Opakování vedlejších účinků

Z claim mechanismu zachovat uživatelský důsledek: procesor má být idempotentní nebo transakční, exactly-once účinky neslibujeme.

[Celá stránka](/home/logik/skola/kkr/kkrunner/docs/persistent-runs.rst)

> Processing ownership uses an OS lock and a persisted attempt token. Observing
> or refreshing a record does not take ownership. A crash can cause a later
> explicit retry to repeat application side effects, so processors should be
> idempotent or use their own transactions. Exactly-once external side effects
> are not guaranteed. The returned Python object is not serialized in the record;
> applications reconstruct historical results from metadata and verified files.

## 6. Životní cyklus a konkrétní akce

R10, R13: stav UNKNOWN, odpojení versus zrušení, attach, načtení historie, retry downloadu a retry parsování.

[Celá stránka](/home/logik/skola/kkr/kkrunner/docs/persistent-runs.rst)

> The normal sequence is ``new -> submitting -> queued/running -> downloading``
> when needed, then ``postprocessing`` when required, then ``finished``.
> ``failed`` and ``cancelled`` are terminal. ``unknown`` means observation is
> uncertain; it is not evidence that the backend stopped.
>
> * ``run.stop_observation()`` detaches this observer without cancelling execution.
> * ``run.cancel()`` requests backend cancellation.
> * ``manager.refresh()`` reconciles backend facts and durable terminal receipts.
> * ``manager.attach(id, process_result=...)`` resumes a nonterminal run.
> * ``manager.resume_download(id)`` reopens a recoverable failed download; attach
>   afterward to resume transfer and completion.
> * ``manager.retry_postprocessing(id, processor)`` retries failed/interrupted
>   processing without submitting the command again.
>
> Do not attach a terminal record. Load its historical application result, or use
> the explicit retry operation for a supported failure phase. A recoverable failed
> download retains its backend outcome receipt until retry settles. Other terminal
> outcomes are acknowledged after the durable state is committed; refresh can
> retry a lost acknowledgement.

## 7. Identita a metadata

R02, R04, R07: jediné ID napříč shardy, řízené aktualizace a včasné hlášení poškozených dat.

[Celá stránka](/home/logik/skola/kkr/kkrunner/docs/persistent-runs.rst)

> ``RunStore.save`` creates a globally unique ID or repeats an identical save.
> Use its atomic update API for subsequent changes. A date shard does not create
> a separate ID namespace. Invalid or duplicate records raise
> :class:`~kkrunner.CorruptRunRecord`; they are not silently ignored.
>
> Metadata is namespaced JSON data, recursively immutable after record creation.
> It cannot change run identity, choose parsers, or execute code. Put the
> application's reconstruction parameters in its namespace and keep executable
> callbacks in the application process.

## 8. Platformy a instalace na cíli

A02: rozlišit čistý Python wheel od skutečné podpory spravovaného backendu; Linux a závislosti cíle.

[Celá stránka](/home/logik/skola/kkr/kkrunner/docs/installation.rst)

> Managed local execution and the direct-process remote service require Linux:
> they use ``/proc``, Unix sockets and process groups for durable supervision.
> Raw local execution uses subprocesses and has no restart guarantee. SSH clients
> need the system ``ssh`` executable. A client does not need ASE or scientific
> parsers even when the application that submitted its runs uses them.

## 9. Ověřování distribučního artefaktu

B01–B04: reprodukovatelný build a test instalovaného wheelu bez zdrojového PYTHONPATH, nikoli staré počty testů.

[Celá stránka](/home/logik/skola/kkr/kkrunner/docs/development.rst)

> The default build first produces an sdist, then builds its wheel. The source
> archive includes documentation sources, tests and Docker fixtures. Build and
> install the wheel in a separate environment; run the tests from the unpacked
> source archive so an editable install cannot hide missing distribution files.
>
> .. code-block:: console
>
>    python -m pip install '.[test]'
>    python -m pytest -q
>    KKRUNNER_RUN_DOCKER_TESTS=1 python -m pytest -q tests/test_ssh_standalone.py tests/test_psij_slurm.py

## 10. Aplikační pravidla patří do ASE2SPRKKR

Do nového aplikačního návodu patří zachování projektových definic, konstrukce
projektového manageru přes aplikační storage policy, výběr vlastního namespace,
obnovení správného readeru, kontrola historických souborů a trvalé pracovní cesty.
Patří sem i cesty v generovaném vstupu relativní k pracovnímu adresáři a správná
báze TAUPATH/TAUDIR. Obecný runner nebude vysvětlovat vědecké registry SPR-KKR.

[Celé navržené znění návodu ASE2SPRKKR](/home/logik/skola/kkr/ase2sprkkr/doc_src/execution.rst)

> Paths written into generated SPR-KKR input are relative to the run's working
> directory, even if the input file itself is elsewhere. An absolute native
> ``TAUDIR`` is mapped relative to the effective ``TAUPATH``. Internal open-file
> paths may remain absolute to preserve ownership and cleanup guarantees.

## Co zůstává pouze v remarks

- Původní seznamy R01–R13/A01–A02/B01–B04, priority a tehdejší verdikty.
- Opravné deníky a chronologie rozhodování; zejména již zrušená kompatibilní API.
- Počty testů, commit hashe, konfigurace konkrétních měření a cesty do `/tmp`.
- Původní duplicitní snapshot architektury v ASE2SPRKKR. Oficiální architektura
  KKRunneru má vlastní stránku; ASE2SPRKKR popisuje aplikační hranici a odkazuje na ni.

## Organizace a publikace

KKRunner: `docs/remarks/`; ASE2SPRKKR: `doc_src/remarks/`. Oba adresáře jsou
vyřazeny ze Sphinx navigace a sestavení. Oficiální stránky používají Sphinx a
Read the Docs theme. KKRunnerův Pages workflow kontroluje HTML i příklady a
publikuje až z `develop`; pull requesty web nemění. Lokální návrh není schválením
publikace. GitHub Pages před schválením a nahráním změn nespouštím.

Schválen je tento výběr a znění odkazovaných návodů s uvedenou úpravou
instalační cesty. Nejde o schválení pravdivosti starých AI protokolů.
