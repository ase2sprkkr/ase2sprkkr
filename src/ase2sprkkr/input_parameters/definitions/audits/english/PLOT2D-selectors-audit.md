# PLOT2DPOT / PLOT2DRHO — SPR-KKR 9.7

Selective input tokens, not new independent solvers. Completed together with shared definition _plot2d_selectors.py. Categories for future joint audit: they work but have a large ALE; selective FULLPOT behaviour and ASA densities are broken/not implemented.

## Dispatch
initvar.f:201–205 sets PLOT2DPR(1) by PLOT2DPOT and flag2 by PLOT2DRHO; both maps internal TASK=PLOT2DPR. GEN (gen.f:1466–1478) selects FPPLOT2DPR for FULLPOT and PLOT2DPR_ASA for ASA.

## Essential Real Behavior
ASA route plot2dpr_asa.f does not contain density array reading or conditional resolution PLOT2DPR flag1/2. Creates potential PLY and geometric line outputs (solid names PLY-V-ASA-... and POT-ASA-LINES.ply). PLOT2DRHO is therefore not the density of the fence to ASA. Definition rejects explicit FULLPOT=false for PLOT2DRHO; actual potential must be FULLPOT, flag does not convert it.

FULLPOT routine fpplot2dpr.f:70 first calls unconditionally PLOT2DPR_FP_PLY. This passes through both potential and density arrays and non-wind outputs according to individual PLOT2DPR flags; to 926 calls STOP_REGULAR with run of test version completed. The following common plane-grid calculation is unattainable in this path. Native selectors therefore do not separate PLY potential/density products as promised by their names. This is not just a input validation problem.

The next STOP in fpplot2dpr.f:671 after writing FULLPOT-pot.pla is in a subsequent unreachable routine; it must not be confused with the first termination. A correction of only this STOP would not suffice. V of this subsequent routine is additionally V<-1 cut to -1 before writing PLY; we do not consider it to be currently available output.

## Geometry and defaults
TASK.AVEC/BVEC/DVEC are three-stage vectors in a.u.; both readers defaultsA=(2,0,0), B=[0,2,0), D=0. B is orthogonalized and maintained at non-zero A.B by its original length. Zero/parallel vectors cause division by zero, definition rejects them as well as non-finite geometry.

NA=50 explicit default corresponds to the reader. NB=50 is explicitly written; if absent NB native is calculated nint(NA*|B|/|A|), which can give0 for very small B. Explicit50 is safer for predictability, but it changes aspect-ratio sampling with different long vectors. NA/NB>=1.

PLOT2DPR_FP_PLY has its own range of fixed sampling parameters: NA/NB are not a guarantee of normal rectangular grid output. Support input parameters does not mean that each of them is used in the script test as in the later plane-grid writer.

GRID0/NE1/ImE0 are explicit startup values; fence task does not exercise an energy loop. Absent FULLPOT for PLOT2DPOT leaves native to recognize potential format; PLOT2DRHO writes true. MPI disabled because GEN guard fence tasks does not permit. DMFT/LDAU fence branches native also rejects; new definition do not offer these MODE parameters.

## Tests and deliveries
7 new cases: both roundtrips, mesh metadata, zero/parallel/infinite vectors, NA0 and rejection ASA density. Full input_parameters set115 passed, two original NumPy warnings. Numerical solver has not been started, outputs are not numerically certified.

Patch0056 independent against baseline267512c4, includes both definitions, common helper and tests. Install all3 definitions files together. Their info explicitly describes native test-writer restrictions.
