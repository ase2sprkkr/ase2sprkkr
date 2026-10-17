# Current integration status

**Updated:** 2026-10-17  
**Branch:** `tasks`  
**File-contract HEAD audited:** `032d932b2e05aa3a540ce232c98522b1fb7065c6`  
**Native source baseline:** `PUB9.7_260420_18h13`

The individual Czech/English source audits and `SPRKKR-ASE2SPRKKR-common-audit.md` preserve the historical source-analysis baselines and historical pytest counts. Statements in those reports that task patches were independent or had not been integrated describe the state at the time of those audits; they are **not** the current repository state.

On branch `tasks`, the audited task definitions and the subsequent file-contract corrections are integrated. The current task-level file-contract audit is authoritative for `input_files`, `cache_files`, and `output_files` roles. In particular:

- required TAU inputs are distinguished from optional reusable TAU caches per task;
- the KKRSPEC family `ARPES`/`AIPES`/`BAND`/`CLPES`/`SPLEED` manifests the native `DATASET_SPEC.out` and polarization-sensitive `.spc` result family;
- `CHI` models `{dataset}.tau` and `TKTK_tables` as reusable cache and manifests the persistent one-pass response outputs;
- `VBPES` models the conditional final-state `{dataset}.tau.dir.F_<photon>` cache and `{dataset}.dos.F_<photon>` output;
- `BSF`/`BSFEK`/`BSFKK` share the exact `{dataset}_BLOCHSF.bsf` and `{dataset}_BLOCHSF_spol.bsf` result manifest;
- `EKREL`, `COMPTON`/`POSANI`, `MEPLOT`/`WFPLOT`, `MECHECK`, `GILBERT`/`STONER`, exchange/current tasks, plot selectors, `FSCAT`/`NUCSSCPL`, `CHIXAS`, and `OPTABS`/`BSDOSIJ` have had their source-defined persistent file contracts integrated and regression-covered.

No persistent-file manifest gap identified by the 2026-10-17 file-contract pass remains on the audited public paths. This does **not** mean that every file opened by native code belongs in a manifest: transient scratch files, disabled/unreachable branches, and files belonging to a known-broken reuse path are deliberately excluded. In particular, the broken `CHILANDAU` `RDTAUFIL` path is not advertised as a reusable cache merely because native code refers to a TAU file.

This integration update does **not** alter the physical/native-code classifications in the source audits and is not a numerical SPR-KKR certification. A task whose native implementation was found broken remains broken even when its persistent result filenames are modeled. Historical test counts in the individual reports remain historical; this document does not claim a new full solver run or a new full pytest run.

## Relevant integration sequence

The current state includes, among others, the file-contract commits for VBPES (`a62b181`), COMPTON/POSANI (`f3a4a7a`), KKRSPEC outputs (`53cfd025`), EKREL (`58511929`), CHI/MECHECK (`d1d3ca7`, `3b384964`), plot outputs (`67f641d`, `0575e130`), GILBERT/STONER (`16bc5e85`), exchange/current outputs (`16eb4e35`), FSCAT/NUCSSCPL (`d73770bf`), CHIXAS (`9fd3b560`), OPTABS/BSDOSIJ (`cf5216e3`), and BSF (`032d932`).
