# VBPES — current implementation audit

**Status:** Valence photoemission via kkrgen, not KKRSPEC; final-state TAU depends on USETAUNN/USETFIN. This describes the integrated `tasks-workflow` implementation against SPR-KKR 9.7 (`PUB9.7_260420_18h13`), not a numerical certification.

## Entry

initvar.f:220–255 recognises VBXPS, VBPES and VBUPS as angularly integrated VBPES; ARPES switches a different branch. The new definition only registers canonical VBPES. Alias names are not registered separately. check_input.f:91–94 requires IREL=3; definition rejects NREL/SREL/SP-SREL. FULLPOT is not prohibited here, but in vbpes.f:224 violently sets TASK.MODE=2.

TASK supportje EPHOT, DEFEP0, ME, MODE, EWORK, GEOMETRY, QVEC, QVECPH, TETPH, PHIPH, NSPIN, ImEFIN, USETAUNN and USETFIN. Section ENERGY additional support is ERANGE which sets EMIN=EF-ERANGE/RY_EV. EPHOT, DEFEP0, EWORK and ERANGE are in eV; ImEFIN and ENERGY.ImE in Ry. TETPH and PHIPH are entered in degrees.

The routine reads directly the initial TAU directory and the number of energy from its metadata (vbpes.f:282–291). The availability and compliance of files cannot be verified only by serialization of the input. USETAUNN uses or calculates TAU end states; USETFIN beats and shuts it down. The calculation of the final states is shifted by EPHOT.

## Defaults: omitted against explicit

| Parameter |Native behaviour|New definition| Serialization default |
| --- | --- | --- | --- |
| TASK.EPHOT | 1253,6 eV | 1253,6 | omitted |
| TASK.DEFEP0 | 0 eV | 0 | omitted |
| TASK.ME | ADA | ADA | omitted |
| TASK.MODE |2; FULLPOT always 2| 2 | omitted |
| TASK.GEOMETRY | 0 | 0 | omitted |
| TASK.ImEFIN | Im(ETAB(NE,1)) |None, delegation per source| omitted |
| TASK.NSPIN |-QHAT after selection of geometry| None | omitted |
| TASK.EWORK |variable EWORKEV is not explicitly initialized before calling; setter gets sentinel 9999| None | omitted |
| USETAUNN/USETFIN |false in absence| false | omitted |
| ENERGY.EMAX |EF at native initialisation| None | omitted |
| ENERGY.EMIN |-0,2 Ry before reading ENERGY or override ERANGE| -0,2 Ry |explicit|
| ENERGY.GRID/NE/ImE |without special task-specific setting GRID/NE; the numbers are then read from TAU| 3 / 300 / 0,01 Ry |explicit|

GRID=3 and NE=300 are the choice of the definition for normal photoemission mesh, they are not a claim of native defaults VBPES. Their writing is explicit. For re-use of existing TAU, the mesh shall correspond to its metadata; the fixed count 300 is not automatic content detection TAU. Better future solution is to get metadata from a specific file. On the contrary, ImEFIN must not be replaced by a fixed forgotten default: it is derived from the last initial energy.

EWORK is only used in this routine to calculate the written EVAC and is not transmitted VBPESINT. Sentinel 9999 is not documented physical default. The definition does not expose it as a known value and does not take over 4,2 eV from a different definition ARPES. If meaningful diagnostic EVAC is required, EWORK must be supplied explicitly by the user. Unverified local initialization EWORKEV is a deficiency native source that should be corrected there.

## Findings

1. **BIS is recognised but does not get over control task.** GEN calls VBPES for VBPES and BIS (gen.f:1445), but vbpes.f:174 requires TASK(3:5)=='PES'. BIS ends here STOP. This definition BIS does not expose; it is an error native of the branch.
2. **GRV is not implemented.** ME_DRIVE has a call for GRV and STOP_MESSAGE (me_drive.f:107–116). The definition allows ADA and NAB which have actual branches of calculation. Neither of them has been tested for numerical accuracy.
3. **EPHOT has a number aliases with a surprising range.** NINT(EPHOTEV)==1/2 overwrites the value to 1253,6/1486,6 eV (vbpes.f:199–200). This is not only an exact value 1 or 2: adjacent energy rounded to these numbers will also be rewritten. The definition of behavior documents and does not overwrite input.
4. **QVEC and embedded geometry are poorly composed.** QVEC is normalized but geometry 0/1/2 then sets a single folder QHAT to 1 without zeroing other or new normalization (vbpes.f:208–250). Without QVEC all components are first reset and the built-in geometry is unambiguous. The definition in this combination warns. Zero and non-finite vectors QVEC, QVECPH, NSPIN rejects the standard due to division.
5. **Custom geometry depends on the completeness of input.** For GEOMETRY outside 0/1/2 is required TETPH at present PHIPH; otherwise QVECPH. The source sets KGEO=3/4 according to the found values. The definition offers 3/4, checks the necessary values and warns when overloading QVECPH angles.
6. **DEFEP0 does not perform a declarable spectrum shift.** In VBPES it is merely transferred and written; not used in the compilation ETABFIN nor transmitted VBPESINT. The same applies to EWORK/EVAC. These parameters cannot be confused with functional energy shifts in KKRSPEC.
7. **POLCONV mixes units of returned angles.** In the vector branch calculates TETPHOT=ACOS(-QHAT(3)) without conversion to degrees but PHIPHOT converts them to degrees (vbpes_aux.f:282–289). The QHAT itself shall be compiled correctly from the vector; the impact on subsequent angular diagnostics requires attention.
8. **MPI postprocessing is not started.** VBPES is not in the list tasks allowed after MPI energy loop in gen.f:1355–1368. The definition has mpi=False.
9. **Some outputs can conflict.** Routine opens solid arpes.dat, arpes_spol.dat and spin.matrix.xps. Names TAU end states contain an integer part EPHOT so that different less than energy can share the same name. The match of the existing file must be checked according to the content.

The correlations are not prohibited in the definition without background: VBPES contains explicit DMFT a branch that sets ORBPOLFIN=NONE for the final state. FULLPOT, correlation and numerical accuracy of the spectra have not been verified during actual run SPR-KKR.

## Current file contract

Required inputs, optional cache and produced outputs follow the integrated task declaration, not a standalone patch. Source: `src/ase2sprkkr/input_parameters/definitions/vbpes.py` and the task registry.

## Verification

Definition tests live in `src/ase2sprkkr/input_parameters/test/`. No new full pytest or native numerical run is claimed.
