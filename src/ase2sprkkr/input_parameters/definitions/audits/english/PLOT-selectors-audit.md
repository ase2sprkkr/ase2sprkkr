# PLOTPOT / PLOTRHO / PLOTSFN — SPR-KKR 9.7

Three selective input options of the same internal calculation branch PLOTPRS, supplemented together. It is not about three independent implementations of graphics. The definitions have a common helper _plot_selectors.py; ZIP contains all four sets for definitions/.

## Native mapping
initvar.f:187–200: PLOTPRS turns on all output flags. If PLOTPRS is not present, PLOTPOT turns on flag1, PLOTRHO flag2, PLOTSFN flag3; any active flag sets internal TASK=PLOTPRS. These selectors can be combined in native input, new individual definitions provide separate canonical options.

GEN patch gen.f:1461–1463 calls FPPLOTPRS. Routine fpplotprs.f:8–11 and the relevant branches provide:
- PLOTPOT: stored VT/BT, for FULLPOT additionally VNST/BNST;
- PLOTRHO: stored RHOCHR/RHOSPN, for FULLPOT additionally RHO2NS;
- PLOTSFN: stored radial coefficients FLMSF and LMISF shape functions for individual meshes, only FULLPOT (r.392). It is not a generator SFN or directly by 3D cell visualization.

It does not perform SCF or DOS. The existing potential/density/SFN data, quality and availability depends on POTFIL. Native outputs are Grace/XMGR, not images created ASE2SPRKKR.

## Defaults and restrictions
TASK and ADSI have an explicit corresponding selectionr. GRID0, NE1, ImE0 are explicit startup values, not claim native fence defaults; energy loop for PLOTPRS is disabled. EMIN/EMAX are not exposed, this branch does not need them.

PLOTSFN explicitly writes SCF.FULLPOT=true and rejects its false. Flag does not change ASA potential to correct FULLPOT file; actual format and SFN data must be present in potential. For ASA native PLOTSFN does not quietly create SFN graphs. PLOTPOT and PLOTRHO leave absent SCF.FULLPOT to recognize file format, work for ASA and FULLPOT.

GEN forbids the PLOTPRS branches with DMFT/LDAU (gen.f:215–223). New definitions do not offer these MODE options; the presence/validity of stored correlation data potential has not been numerically certified. MPI is disabled, native GEN MPI guard PLOTPRS does not permit. If some types or meshes lack the corresponding data, the input definition alone will not detect this.

Category for future joint audit: have fully meaningful as selective diagnostic output; PLOTSFN conditional on actual FULLPOT potential. It's not the missing physical solver.

## Verification
4 new regression cases: roundtrip and pickle each definition, serial/executable metadata, explicitly entered FULLPOT for PLOTSFN and rejection false. Full input_parameters set 112 passed, two original NumPy warnings. Numerical native run not started.

Patch0055 is independent from baseline267512c4; does not require earlier PLOTPRS patch. All four definitions files must be installed together.
