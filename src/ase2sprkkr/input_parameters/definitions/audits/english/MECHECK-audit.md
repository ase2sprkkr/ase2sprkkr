# MECHECK — current task audit

`mecheck.py` exposes the serial `kkrgen` matrix-element diagnostic `ME_TEST_11` selected by `initvar.f:162–167` and called from `gen.f:1458–1459`. Limitations of the compiled SPR-KKR 9.7 routine remain even when the Python input definition accepts the parameters.

## Interface

MECHECK is a diagnosis matrix elements transmitted `gen.f:1455–1459` to `ME_TEST_11`. `initvar.f:162–167` sets ITEST=11 and requires KKRGEN. The definition uses `kkrgen`, MPI=False: the routine does not have a division of its work between the wounds and writes shared file names. Numerical start-up or behaviour MPI has not been verified.

| Parameter TASK |Default when missing/behavior|Definitions|
|---|---|---|
| ITSEL |1; Atomic type index, not site|positive integer, 1 is omitted from the input|
| PRINT |preserves CONTROL.PRINT|dynamic default, same value is omitted from the input; here actually PRINT, not MEPLOT.IPRINT|
| INITIALSTATE | BAND |BAND or CORE; default is omitted from the input|
| CL |2P for CORE| shell n=1..6, S/P/D/F; validace l<n; default is omitted from the input |
| ERYDI |BAND: effectively {0.5,0.01} Ry|Exactly two real components, default is omitted from the input|
| ERYDF |BAND: effectively {0.8,0.01} Ry|Exactly two real components, default is omitted from the input|
| NONMAG | false |TASK flag suppresses BT of the selected type; different from CONTROL.NONMAG|
| NMEFORM |Reader starts at 3 but the effective value is always 2|Reader default 3 is omitted from the input; explicit setting will indicate that overwrites|

Source: `me_test_11.f:131–177, 297–346`. BAND evaluates only one pair of complex energies; CORE calculates initial state from CL and final states from the first grid ENERGY. ERYDF is readable for CORE but overwrites ETAB; ERYDI is not loaded in CORE. The explicit ignored choice of definition is highlighted.

If the absolute value BAND of the energy is less than 1e-6, the source shall replace it with the above default. Therefore, it is not possible to demand exactly zero with this input. The definition draws attention and does not overwrite the input. Fortran reader also accepts one real component and copies it to another (`scan_input.f:698–707`): scalar here does not mean complex energy with zero imaginary part. The definition requires an explicit pair not to be hidden.

## Restrictions

- `check_input.f:99–100` requires IREL=3. NREL/SREL/SP-SREL definition rejects.
- `gen.f:211–222` explicitly rejects DMFT and LDAU and for MECHECK. The same restriction is in the definition.
- FULLPOT is generally not prohibited for MECHECK in the source; the frequency calls FPSSITE and has corresponding branches. The definition therefore does not invent a ban. Its numerical accuracy is not verified here; the mode is derived from potential.
- SITES.NL remains None and is omitted from the input: cutoff determines SPR-KKR/potential. U explicit uniform cutoff is checked for availability CORE orbital. The non-uniform cutoff and the availability of a particular atomic-type shell without potential cannot be verified.

## Problems in resources SPR-KKR

1. **NMEFORM has no effect.** After loading and possibly FULLPOT the limit follows `NMEFORM=3` and at a fixed `CMP_NAB_VS_NAB=.FALSE.` additional overwrite to 2 (`me_test_11.f:96, 172–177`). Running NAB and ADA; GRV cannot be switched on by input.
2. **The second check described in the comment does not run.** Compilation constant `K_CALC_ME=2` sets NTEST=1 (`:49–50, 370–374`). Replacement J using Z and subsequent checking of irregular members are unavailable; this is not a missing input parameter.
3. **Summary statistics are not written.** Loop `DO IMEFORM = 2, - NMEFORM` (`:636`) with default step +1 has no iteration at NMEFORM=2. Probable typo in the top between; partial comparison before it runs.
4. **Side writing to unit 601.*** NAB to branch CMATSTRUCT with unit 601 (`:478–480`) for which this routine does not contain OPEN. The default file `fort.601` can be created according to Fortran runtime. It's a suspicious remainder of the diagnostic, not the input parameter; runtime impact hasn't been verified.

## Defaults ENERGY

ASE2SPRKKR here explicitly writes GRID=3, NE=1, ImE=0.01 and EMIN=-0.2 Ry. These are deliberate choices of simple diagnostics, not claims of native defaults. For BAND the routine is ignored when selecting both state energy. For CORE actually define the final energy path; for spectrum it is necessary to increase NE and set appropriate limits. NE=1 is suitable for basic control, for CORE scan would be higher default is more practical; change is a design of user behaviour, not repair of parser. EMAX=None is omitted from the input and delegates actual Fermi energy SPR-KKR. Couple EMINEV/EMAXEV can be used for a specific window.

## Native output contract

Compiled `K_CALC_ME=2` makes the named `DATSET_check_*.dat` diagnostic files unreachable: native `me_test_11.f:403–411` opens them only if `K_CALC_ME>2`. For a CORE energy scan with more than one final-state energy and a nonzero initial matrix element, `me_test_11.f:704–742` instead opens zero or more `MPZAZB_{-,z,+}_K-I-J.dat` files. These names **do not contain DATASET or ADSI**. Current `mecheck.py::output_files()` correctly declares wildcard families only for `INITIALSTATE=CORE` and `ENERGY.NE[0]>1`; the exact count is data-dependent. The routine may also write to unnamed Fortran diagnostic units, but their runtime names are not stable canonical output files.

## Verification scope

Input serialization and validation have dedicated regression tests; no fresh test-suite run or numerical SPR-KKR calculation is claimed here.
