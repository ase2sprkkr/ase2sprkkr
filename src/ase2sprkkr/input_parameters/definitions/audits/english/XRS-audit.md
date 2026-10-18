# XRS © experimental definitions and audit

SPR-KKR 9.7 from sources supplied versus ASE2SPRKKR refactor after fetchi, commit 267512c4b8ccbbfd37af7ede898febc066b7e3ee. XRS has no definition yet. Single patch 0026-xrs.patch adds xrs.py and test_xrs.py; does not require XAS/XES patche.

## Real interface

XRS calculates complex tensor resonance X-ray dispersion amplitudes in a common routine XRAYSPEC. Unlike RXES it has its own factors, resonance denominators and output SIG. Input TASK is XRS. CHECK_INPUT requires FREL. INITVAR activates IXRAY, reading TAU and energy loop.

For XRS three energy integrations are used: the first with TAU and the reading of the single-site part, the second with the single-site part on the first network's real energy and the third with EPATH GRID7 from ENERGY.EMAX to TASK.EMAX3 with TASK.NE3 points (xrayspec.f:684–702,1213–1220). TSELECT=0 includes both members; 1 resets non-resonance WTERM2, 2 resonance WTERM1 (413–419).

TASK contains IT, CL, ME, POLOBS, DIPOL, GFCONVRH, WRREGIRR, non-dipole flags, three Euler frames and TAUCORE. The widths TAUCORE are in eV, not more than two elements for core kappa branches; the common reader repeats the last value. Values below1e-6 take WCOREHOLE. Unlike XAS/XES, they have proven use in resonance denominators (1259,1264).

OUTPUT is not offered: XRS rewrites it to SIGMA regardless of input. MESPLIT is converted to true. DYNAMICAL is not read for XRS. These options would not represent functional switches task.

## Defaults and writing

| Option |Native missed effect|Definitions|
| --- | --- | --- |
| IT / CL | 1 / 2P |Same values, omitted|
| ME / POLOBS / DIPOL | ADA / S / 0 |Same values, omitted|
| NE3 / EMAX3 / TSELECT | 40 / 8 Ry / 0 |Same values, omitted|
| TAUCORE |[0,0], delegation to WCOREHOLE| None, omitted |
|Frame angles|0/0/0, absent CHANGE_FRAME=false| None, omitted |
| GFCONVRH / WRREGIRR |false by Keyword absence| false, omitted |
| ENERGY.GRID / NE |XRS has no own default branch; NE begins0|explicitly6/180|
| ENERGY.EMIN |XRS then sets to EF|omitted, section has no lower limit|
| ENERGY.EMAX |general overwrite of the sentinel to EF|explicitly4 Ry|
| ENERGY.ImE | 0,01 Ry |explicitly0,01 Ry|

GRID6/180 and upper limit4 Ry represent the explicitly selected range of vacant stocks similar to native XAS scan. XRS alone does not have native EMAX4; its general EF would lead to the first interval EF..EF, therefore it is not more suitable for this scan. Deleted EMIN takes over the actual EF and does not represent an ineffective absolute default. The convergence of the first and third intervals has not been numerically verified. ImE corresponds to native value.

Validator controls FREL, physical core shell n<=5/l<n, positive IT/NE3, NE>=2, one energy path, finish angles and widths, non-negative TAUCORE and EMAX3>ENERGY.EMAX at known absolute between. Upper limit IT, existence shells and conformity ENERGY.NE with the first TAU contour require actual potential and TAU directors.

## Serious errors and restrictions

1. **SIG is initialized only partially:** allocation SIG(NEMAX+50,NCSTMAX,3,3) at line 636 but CINIT has a length NEMAX*NCSTMAX*3*3 at line 639. It is then added over NESIG=50+NE. Due to the column saving, these are different areas of core/polarization channels, not only the last50 energy. The result may contain uninitialized data; the correction requires Fortran.
2. **The optional network is not monotonous:*** before the first energy point is generated50 points by the expression (NESIG-IE-1)^3. IE=49 already corresponds to the first energy and IE=50 lies slightly above it, then reconnects the first energy TAU. Metadata of both number and axis sorting must respect the result reader.
3. **MPI does not occur spectrum:** GEN:1355–1368 for XRS ends before calling XRAYSPEC; definition has mpi=False.
4. **The .rat filename and header retain the X-ray absorption conventions:** XRS can collide with XAS/XES for the same DATASET, type, and shell. Output XRS writes for each energy and core state9 complex components SIG. It cannot be interpreted as3 absorption polarisation intensity.
5. **The output number of points is not simple ENERGY.NE:** NESIG starts at50 and connects NE; at the same time, the later third integration changes NE to NE3. The reader shall not derive the number of rows from the current NE or its original value without knowledge of the structure.
6. **Non-dipole flags do not guarantee non-dipole XRS tensor:** their field is filled but owns XRS SIG the sum uses MZAZB/MZBZA and their conjugate variants. Based flags correspond to a common reader; the effect on the full dispersal tensor requires additional physical audit. ME can also not be used as verified independent implementation of three Gauges.
7. **Reserve for third integration:** if NE3>NEMAX, XRAYSPEC reallocates ETAB/WETAB/ETAB0. The next previously allocated fields will not automatically grow. Access to these fields in the third branch must be verified before using NE3 above the original capacity; this patch allocation errors does not fix.

Contrary to FORCETEN it is not a setholder off routine, therefore the input was supplemented. Due to partial initialization SIG, the planned joint audit includes the broken tasks for the repair of Fortran. No numerical result is certified by this definition.

## Tests

16 new cases: serialization/load/pickle, omitted native defaults and explicit EMAX, resonance parameters and TAUCORE, relativity, core shells, split contour, IT/NE3/TSELECT, finish and order of the third limit. Full input_parameter set: **124 passed**. SPR-KKR not started. Patch passed git diff --check and re-application control.
