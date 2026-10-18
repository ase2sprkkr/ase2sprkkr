# PLOTJALF / PLOTJORB — SPR-KKR 9.7

Category for joint audit: **broken**, not generally disabled placeholder. This is the postprocessing of options JALF/JORB, not the conversion of currents.

## Intended entry
CONTROL.PLOTJALF and CONTROL.PLOTJORB are actual flags (initvar.f:122–123), not separate TASK keywords. Each set internal TASK=PLOT2DPR (lines 205–206). The geometry is read from TASK.AVEC/BVEC/DVEC/NA/NB. In the expected current branch fpplot2dpr.f:194–231 is loaded DATSET0_JALF.dat or DATSET0_JORB.dat; if both are set simultaneously, JALF has priority. It requires corresponding types, radial meshes and cutoff, otherwise STOP data file not consistent. The input data should come from an earlier JALF/JORB calculation.

## Why the definition cannot ensure functionality
GEN Dispatch (gen.f:1466–1478) for internal PLOT2DPR select FPPLOT2DPR to FULLPOT, PLOT2DPR_ASA to ASA.

FPPLOT2DPR calls unconditionally PLOT2DPR_FP_PLY (fpplot2dpr.f:70) at first. This routine does not read PLOT_JALF/PLOT_JORB or current file. It works with potential/density and ends at the end STOP_REGULAR with run of test version completed (plot2dpr_fp_ply.f:926). Return to the subsequent current branch is therefore not started in normal FULLPOT even if PLOT2DPR(1/2) is not switched on. Conditional RETURN to 103 does not help here: requires NOT FULLPOT.

ASA routine PLOT2DPR_ASA current flags or data does not read, draws potential. Therefore, the change FULLPOT flag is not bypassing the problem.

The result may be a different plot output or an earlier error, but never this documented path required by rendering saved JALF/JORB field. The presence of the current file will not repair the dispatch itself.

## Delivery and other work
Not creating input definition or Python test copying an unattainable native branch. First it is necessary to separate the test PLYwriter from the current branch in FPPLOT2DPR and verify the numerical postprocessing of the reference JALF/JORB file. Only exposure CONTROL flags in ASE2SPRKKR would not be functional support. Native Fortran was not changed; numerical run did not take place.
