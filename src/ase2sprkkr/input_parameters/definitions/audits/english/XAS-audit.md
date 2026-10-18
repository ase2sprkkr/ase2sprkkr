# XAS audit of existing definitions and corrections

Compared with the supplied SPR-KKR 9.7 and the fetched ASE2SPRKKR refactor, commit 267512c4b8ccbbfd37af7ede898febc066b7e3ee. Patch 0025-xas.patch regulates existing xas.py and adds test_xas.py. He's single, he doesn't need a previous patch.

## Imported errors of the original definition

1. **MPI=true does not complete the spectrum.** GEN:1355–1368 after MPI energy loop stops all tasks outside the explicit list; XAS or XES are not in it. XRAYSPEC is called to 1378–1385. Patch uses mpi=False. Native MPI guard is a source proof, not a numerical test.
2. **Default EMAX was not a real parameter.** ENERGY(defaults={"EMAX":4.0,...}) EMAX does not create: common ENERGY reads the limit from the emax argument, not from this dictionary. The original input could rely on native XAS default4 Ry, but the definition did not allow standard setting ENERGY.EMAX. Patch uses emax=(4.0,...,None), so it can be adjusted and4 Ry is explicitly written out.
3. **MECHECK is not a local control XAS.** INITVAR:162–167 recognizes MECHECK as a separate task/test11 and sets other workflow flags. Later detection XAS is not proof of safe combined mode. XRAYSPEC does not read this option. Patch removes it from XAS TASK; separate definition MECHECK is an earlier separate patch.
4. **CL as free str did not watch the physical shell nor the table boundary.** XRAYSPEC has SHELL(5) s K,L,M,N,O and directly indexes NCXRAY at line 491. n=6 would cross the field. Newly authorised n=1..5, s/p/d/f and l<n. Existence shells in a specific potential is verified only outside this validator.
5. **IT as an int allowed0 even negative numbers.** Patch requires a positive index type. The upper limit requires potential content.
6. **Frame labels confused the rotation of the frame as one light vector.** INIT_MROT_FRAME reads three Euler angles and forms GETMROT(-gamma,-theta,-phi). This is not a simple parameter of one electric field; XRAYSPEC counts three polarization channels. Labels have been corrected and missing FRAMEGAM supplemented.

## Parameters and validation

Completed ME, POLOBS, DIPOL, GFCONVRH, WRREGIRR, NONDIPE2 and NONDIPM1 by XRAYSPEC:295–329,442–465. POLOBS=S vypisuje (+,-,z), C vypisuje (x,y,z). ME=NAB overwrites DIPOL to 1. If both non-dipole flags are present, NONDIPE2 takes precedence. The definition draws attention to both combinations.

The validation rejects other MODE than FREL, invalid core shells, non-quantifiable framework angles, inexplicable/single point NE and more energy routes. CHECK_INPUT:43–47 directly requires IREL=3. XRAYSPEC:376–380 requires compliance ENERGY.NE with the number of points of the first TAU contour; content TAU and the potential of the module does not read, so it does not issue this match as verified.

Not added MESPLIT (after reading is K_ME_SPLIT always true), TAUCORE as supposedly effective management XAS breaking or DYNAMICAL without auditing its separate routine. TAUCORE with a common reader reads but its effective application is in branches XMO/XRS; standard call ME_ALF_ALF sets CECOR=0. DYNAMICAL for RXAS/XMO actually switches to XRAYSPECDYN, which is not this audited static workflow.

## Defaults: written to the input and omitted

| Option | Native XAS |Corrected definition|
| --- | --- | --- |
| ENERGY.GRID | 6 |6 explicitly, unchanged|
| ENERGY.NE | 180 |180 explicitly, unchanged|
| ENERGY.ImE | 0,01 Ry |0,01 explicitly, unchanged|
| ENERGY.EMIN |After a general reader XAS default EF|-0,2 Ry explicitly, left|
| ENERGY.EMAX | 4 Ry |4 Ry explicitly, newly actual parameter|
| TASK.IT / CL | 1 / 2P |1 / 2P explicitly, left|
| OUTPUT |BARN, output in Mb|None omitted; accepts BARN, MBARN, SIGMA|
| POLOBS | S | S omitted |
| ME / DIPOL | ADA / 0 | ADA / 0 omitted |
| GFCONVRH |false in absence of keyword| false omitted |
| WRREGIRR |Effectively false: SECTION_FIND_KEYWORD overwrites DATA true| false omitted |
|Frame angles|0/0/0, without keyword CHANGE_FRAME=false|Same zeros omitted|
| MODE | FREL | None omitted |

EMIN=-0,2 Ry is not erroneously omitted default: writes se. Native EF is usually more natural for the standard absorption interval and avoids parts below EF. I don't change it automatically because existing ASE2SPRKKR behavior can serve a deliberate scan of the wider interval. The specific suitability depends on EF and the required spectrum; convergence has not been tested.

Note: in INIT_MOD_ENERGY:234–249 is XAS in absence EMF returns to EMIN=EF, so the absolute EMIN itself may not always remain effective. Relative pairs EMINEV/EMAXEV and its subsequent processing require inspection of this branch. Patch does not introduce a fictitious guarantee that Python default always determines the actual start native scan.

The output SIGMA is in10^15/s, the other OUTPUT uses Mb. Original DefKeyword had no explicit base default; missed OUTPUT correctly left native BARN. New Keyword keeps this missed and adds directly native name BARN next to compatible MBARN.

## Other native restrictions

- MEFORM is loaded, NAB affects DIPOL and GRV prepares the derivative. The standard calculation still calls ME_ALF_ALF without giving MEFORM. It cannot be argued that three accepted names mean three independent, numerically verified Gauge implementations.
- Non-dipole call is only for .NOT.FULLPOT and DIPOL=0. FULLPOT does not therefore need to evaluate the same contributions with the same flags.
- The calculation selects a single type IT (NTXRAY=1). It does not issue a list of types.
- The number of points, mesh and type of potential must be consistent with TAU directors. Only successful serialization input does not test this condition.
- XAS and XES with the same DATASET, type and skin may use the same .rat name. It is necessary to distinguish outputs at work with both tasks.

## Tests

15 new cases: defaults and explicitly written to the input EMAX, save/from_file/pickle, own EMAX and frame/OUTPUT/POLOBS, elimination MECHECK, core shell limits, IT and DIPOL, relativity and split contour, non-finite angle, NAB override and non-dipole flags priority. Full input_parameter set: **123 passed**. Numerical SPR-KKR and the accuracy of the .rat spectrum have not been tested. Patch passed git diff --check and re-application control.
