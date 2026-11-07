# Execution boundary and public contract

KKRunner owns ExecutionSpec/ExecutionOutcome, TaskFile, ExecutorEndpoint,
RunState/FailurePhase, Execution/WorkingDirectory, RunRecord/RunStore,
FileEntry/FileMapping/FileManifest, FileState/FileVersion/FileFingerprint,
ExecutionResources and ResourceOwnedRunner. Public exports live in `kkrunner`.
ProjectRepository and the public ExecutionManager organize tracked runs in
projects; `kkrunner.manager.ExecutionManager` also supports an explicit RunStore.

ASE2SPRKKR owns input preparation, executable selection, task definitions,
template/callback resolution for TaskFile declarations, parsers, TaskResult and
the `ase2sprkkr` application-metadata namespace. It depends on KKRunner.
KKRunner depends only on platformdirs, plus optional psij-python.

TaskFile is an immutable declaration. Its `type` and `info` are opaque strings;
KKRunner does not resolve application templates or consult format registries.
FileEntry is a resolved persistent location plus transfer/publication state.
FileManifest owns the authoritative declarations, cache baselines and captured
versions. Read/write conflicts derive exclusively from the manifest; task names
have no locking meaning. Clients must declare every read/write/cache footprint.

RunRecord owns persistent run identity, backend identity, manifest, namespaced
JSON metadata, `processing_required` and the latest processing-attempt ID.
Metadata is recursively immutable and cannot change through RunStore.update.
RunRecord schema_version is 1 throughout pre-release development, even when
its structure changes. Only the current shape is accepted. Development records
are not migrated. Old ASE2SPRKKR schema-3 records are intentionally unsupported;
copying their JSON into the new store raises CorruptRunRecord. User configuration
and Project generation/versioning are separate contracts.

Lifecycle: new -> submitting -> queued/running -> downloading (if required) ->
postprocessing (if required) -> finished. Any phase can fail or be cancelled as
supported by the backend. UNKNOWN is recoverable uncertainty, never completion.
Download and processing failures can be explicitly reopened without resubmission.
Without a required callback, successful execution remains POSTPROCESSING.
The callback return object is not persisted; historical applications reconstruct
it from the recorded application metadata and verified files.

Postprocessing claims a per-run filesystem lock, captures/verifies artifacts,
publishes the attempt, calls the processor and commits success/failure before
releasing ownership. Other clients may monitor but cannot process concurrently.
OS lock release after process death makes an interrupted POSTPROCESSING attempt
recoverable. A live owner is never displaced based on a timeout. All clients of
a store must share the same reliable file-locking namespace. Callbacks may be
repeated after a crash and must tolerate that (idempotence or transactional side
effects); exactly-once side effects are not promised.

Backend facts and durable submission/outcome receipts are not a second run
record. RunRecord remains the authority for end-to-end completion. Remote ACK
follows its durable commit. JSON writes use atomic replacement, file and directory
fsync on POSIX; portable Windows directory durability is weaker as documented
in _durability.py. Invalid records raise errors, including during enumeration.

Before submission the client/prepared calculation owns its open resources.
Successful create_run transfers them to ResourceOwnedRunner. Creation/submission
failure closes acquired resources, completion closes the owned resources, and
cancellation before submission closes them without launching. Borrowed handles
retain their existing ownership semantics. Managed execution uses file paths on
the execution host; the supervisor owns process lifetime beyond client detach.
Temporary directories must remain available for backend execution and processing;
clients must not delete them on mere loss of observation.

The CLI is `kkrunner remote`. The original transport framing/control protocol is
preserved and still rejects incompatible protocol versions. Worker module paths,
default command, runtime namespace and PSI/J identity markers now belong to
KKRunner. A remote host needs only KKRunner and the desired executable (plus
PSI/J for scheduler use). `ase2sprkkr remote` delegates only the end-user status/follow/cancel subset.
`kkrunner-remote` is an alternate entry point for the same protocol implementation.

Migration is split into extraction, metadata/API integration, and recovery changes.
No duplicate runner implementations or result-processor registry is introduced.

Tracked ASE2SPRKKR runs reject automatically deleted working directories
(`directory=False`) and temporary stderr (`error_file=False`); callers select
persistent paths to make restart recovery and history possible. Ordinary
untracked temporary runs keep their previous cleanup behavior.

TaskFile.effective_info is no longer part of the data type; ASE2SPRKKR clients
can use result_file_definition_registry.info_for_task_file(declaration).
Application-manager load_result/attach_sprkkr are replaced by the existing
ASE2SPRKKR outputs functions load_project_result/attach_project_run.
