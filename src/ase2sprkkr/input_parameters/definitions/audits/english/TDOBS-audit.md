# TDOBS — SPR-KKR 9.7 / ASE2SPRKKR

New definition TDOBS for kkrgen, serial. Category: works, but has major caveats © is not a separate new analysis compared to generic DOS/observables.

## Testified behaviour
initvar.f:523–530 detects TDOBS, sets RDTAUDIR=true, WRTAUDIR=false and turns on initialization/energy loop. kkrgen.f:932–961 checks tau cache: with available cache initially shuts down RUNELOOP/INITELOOP, without it sets recuperation and writing tau. Then kkrgen.f:971–974 for TDOBS always turns both loops back on. GEN in energy loop calling CALCDOS (gen.f:869–876) and CALCDOS_WRITE (r.1151). No other dispatch branch TDOBS was found in the sources.

CALC_VECOBS is a different internal switch, initialized false in mod_calcmode.f:116; parser TDOBS does not set it. The name TDOBS is therefore not represented as a guarantee of vector or time-dependent analysis. General DOBS areays themselves do not support such analysis.

GEN MPI guard (gen.f:1360–1371) TDOBS does not permit and terminate STOP_REGULAR; mpi=False is required. Support general DOS energy loop is not numerical verification of all cases CPA/FULLPOT or re-use of cache.

## Inputs and defaults
TASK=TDOBS and CONTROL.ADSI=TDOBS are recorded. TAU.TAUDIR/TAUPATH are optional without Python default: native uses DATASET.tau.dir and current directory (kkrgen.f:266–273). Cache must correspond to potential, energy and configuration; definition does not control its content.

GRID3, NE300, ImE.01, EMIN-.2, EMAX1 Ry are explicit spectral defaults taken from audited DOS definition. They do not have to correspond to generic native initial defaults but are not silently omitted. EMAX1 is not EF sentinel and includes vacant states; a narrower window may be more suitable according to the objective. We do not exhibit RDTAUDIR as an untrue overwritten option: its true is native behavior task. The missing tau does not automatically cause the end; native selects the recomputation fallback.

Validation NE>=1, GRID2/3 at least2 points, final bounds, ImE>=0 and warning at GRID3/ImE0. Added MODE, CPA and STRCONST as well as for DOS. The definition is independent of the previous DOS patchi; it does not inherit any non-application correction.

## Verification
5 new cases: roundtrip TASK/ADSI and tau paths, MPI/executable metadata and rejection of invalid mesh/bounds. Full input_parameters set 113 passed, two original NumPy warnings. Numerical SPR-KKR not started, cache reuse has not been numerically tested.
