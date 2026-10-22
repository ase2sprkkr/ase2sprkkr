# T1 audit and new definition InputParameters

Background: attached SPR-KKR 9.7 (PUB9.7_260420_18h13) and current ASE2SPRKKR `refactor`, commit `267512c4b8ccbbfd37af7ede898febc066b7e3ee`, verified fetch. T1 does not have a domain module here. The separate patch adds `definitions/t1.py` and `test/test_t1.py`; it is not dependent on previous patches.

## Entry and modes

`initvar.f:331–334` recognizes T1 and turns off RUNELOOP. GEN at T1 rejects FULLPOT and calls NONRELT1 for IREL<=1, otherwise RELT1 (`gen.f:1485–1496`). Routines do not read task-specific parameters. Therefore, I do not create ineffective options ITSEL, temperature, isotope or magnetic field. Nuclear data come from internal tables TABGAMMAN; other tables in non-functional NONRELT1 are not user input parameters.

The definition uses CONTROL, TAU, MODE, STRCONST, ENERGY, TASK, SITES; executable=kkrgen and mpi=False. Postprocessing RELT1 is not distributed between wounds and uses local TAUT. Numeric or MPI startup not verified.

- FREL/None: input for RELT1; default mode.
- NREL and SREL: NONRELT1 contains `IF (NT.GT.0) STOP_MESSAGE('not available')` (`nonrelt1.f:83–84`). The definition rejects these schemes because they are not available for the actual system NT>0.
- SP-SREL: GEN sends it to RELT1 but the routine works with relativistic kappa channels. The definition draws attention to an unverified combination, does not invent explicit native ban.
- DMFT/LDAU: prohibits GEN (`gen.f:211–222`); definitions also.
- FULLPOT: zakazuje check_input and GEN. The definition documents the prohibition; the input alone cannot reliably determine FULLPOT mode without loading the potential.

## Energy and defaults

`init_mod_energy.f:265–275` for T1 sets hard GRID=3, NETAB(1)=1, EMIN=EMAX=EF and ETAB(1,1)=EF+i*EIMAG. EIMAG begins for this task to 0 and can be rewritten ENERGY.ImE. In addition, RELT1 creates CEFERMI=DCMPLX(EFERMI,EIMAG) for SSITE (`relt1.f:159–171`). The energy window and spectral scan are not effective parameters here.

| Option |Definitions| Relationship to SPR-KKR |
|---|---|---|
| ENERGY.GRID |explicit 3|is converted to 3|
| ENERGY.NE |explicit 1|is converted to 1; non-zero initialization required before reset|
| ENERGY.ImE | 0 |equal native default task; user non-zero ImE is actually read|
| ENERGY.EMIN/EMAX |not in the definition|T1 is always overwritten to EF; there is no reason to add misleading interval|
| SITES.NL | None, is omitted from the input |cutoff determines potential; explicit NL>5 rejected due to capacity CHK|
| CONTROL.PRINT |explicit 1|deliberate difference from general 0 due to initialization error NN/MM|
| CONTROL.RUNELOOP |explicit flag true|overwrites native loop off to create TAUT|

GRID/NE different from {3}/{1} definition only warns as overwritten. ENERGY.NE must still be positive: INIT_MOD_ENERGY allocates ETAB and controls NEMAX>0 before task-specific T1 reset.

**PRINT=1 and RUNELOOP must be published.** This is not omitted defaults which should correspond to native: these are explicit options ASE2SPRKKR passing by two initialization problems listed below. V in this case native defaults are not better. Definition at PRINT=0 or off RUNELOOP rejects saving. Input for use of already prepared TAU data with other loop control does not solve this definition.

## Serious errors and limitations SPR-KKR

1. **Not initialized NN/MM at PRINT=0.*** RELT1 sets NN=2*NL² and MM=NKMMAX only inside `IF (IPRINT.GE.1)` (`relt1.f:150–154`), but later always transmits them to CMATSTRUCT (`:187`). General PRINT=0 may therefore lead to incorrect access to memory in the diagnostic extract. Explicit PRINT=1 in the definition ensures initialization; RELT1 then sets it itself to 0.
2. **T1 does not calculate normally TAUT which it then uses.** INITVAR turns off RUNELOOP for T1 (`:331–334`). RELT1 counts SSITE again, but it solves single-site TSST, not TAUT; uses `W2(:,:)=TAUT(:,:,IT)` (`relt1.f:206`). TAUT is obtained in an energy loop GEN which is `IF (RUNELOOP)` (`gen.f:534`). For separate T1 run without prepared TAU matrix it is a faulty initialization path. Explicit CONTROL.RUNELOOP activates later native override (`initvar.f:689–693`) and allows a loop at one point EF. This procedure is verified by reading resources, not numerically during. When manually added RDTAU and existing TAU files, the main program can change flags again; this mode has not been checked.
3. **The negative kappa as the upper limit of the loop.** V the sum of the relaxation contributions is calculated signed KAP0/KAP1/... and loops are used `DO M0=1,2*KAP0` and similar to M1/M2/M3 (`relt1.f:215–245`). Negative kappa channels have no iteration; for example with1/2 withkappa=-1 will be omitted. The number of m_j states is related to the absolute value of kappa. It is not just missing posts: shifted LAM indices then do not correspond to the whole kappa base. The very definition of input does not correct this numerical error.
4. **Restricted label capacity.** CHK has 10 items (`relt1.f:39,54–57`), but extracts use NKQ=2*NLQ-1 (`:338–347`). For NLQ>=6 accesses an array out of bounds. The definition rejects explicit max(SITES.NL)>5; at omittedNL, the cutoff must be checked in potential. I do not declare the implicit context safe without loading it.
5. **NONRELT1 is intentionally blocked.** Its extensive code for STOP is not proof of an available nonrelativistic calculation. The definition does not represent these branches as supported.

Before using relaxation speeds for physical conclusions it is particularly necessary to correct signed kappa loops and compare numerical results. Patch adds input interface and initialization protection, not validated solver T1.

## Output and tests

RELT1 writes tables and diagnostic matrix into standard output. The specialized parser of these tables is not part of the patch; the generic output reader is used.

10 new tests: Discovery, explicit PRINT/RUNELOOP, file/text round trip, pickle/copy, non-enabled energy limits, unavailable NREL/SREL, forbidden correlation, PRINT=0 and off loop, transcription GRID/NE, kappa label capacity and efficient ImE.

Whole set `src/ase2sprkkr/input_parameters/test`: **118 passed** to given default revision plus this separate patch (without previous patches). `git diff --check` and the inverted application of patch has also passed. Numerical run SPR-KKR not performed; tests verify input definition.
