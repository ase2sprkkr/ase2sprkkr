# PLOT2DPR audit and new definition InputParameters

Background: attached SPR-KKR 9.7 (PUB9.7_260420_18h13) and current ASE2SPRKKR `refactor`, commit `267512c4b8ccbbfd37af7ede898febc066b7e3ee`, verified fetch. PLOT2DPR does not have a definition here. Patch adds `definitions/plot2dpr.py` and `test/test_plot2dpr.py` separately; previous patches are not necessary.

## Input and variants

`initvar.f:199–205` recognises PLOT2DPOT, PLOT2DRHO and PLOT2DPR. The last keyword will turn on both groups. GEN calls in FULLPOT FPPLOT2DPR, in ASA PLOT2DPR_ASA (`gen.f:1465–1479`).

Creating: `InputParameters.create('PLOT2DPR')`. The selective variant can be set over `p.TASK.TASK = 'PLOT2DPOT'` or `'PLOT2DRHO'`. This patch does not register separate modules for `create('PLOT2DPOT')` and `create('PLOT2DRHO')`; support is their keyword in the definition PLOT2DPR and read input over `from_file`.

| Parameter TASK |Default when missed|Definitions|
|---|---|---|
| AVEC | {2,0,0} |Three real components, default is omitted from the input|
| BVEC | {0,2,0} |Three real components, default is omitted from the input|
| DVEC | {0,0,0} |Three real components, default is omitted from the input|
|NA| 50 |positive number of intervals, default is omitted from the input|
| NB |NINT(NA*length(BVEC)/length(AVEC))|None, is omitted from the input; translated SPR-KKR|

Vectors are cartesian in units ALAT, **not broken coordinates based on grid vectors**. AVEC and BVEC form rectangular window edges; DVEC moves it from a common beginning. In the older regular FULLPOT branch NA/NB indicate intervals and the grid has (NA+1)*[NB+1) vertices. The actual dialed PLY routines also have their own radial/angle sampling parameters firmly in the code: NA/NB cannot therefore be used as a full density check resulting PLY triangulation.

Source defaults: `plot2dpr_asa.f:53–54,64–85`, `plot2dpr_fp_ply.f:76–77,126–149`, `fpplot2dpr.f:63,102–125`.

## Geometric validation

Both branches adjust the non-zero scalar product BVEC by projection perpendicular to AVEC and maintain its original length. The definition highlights this change, the original input does not calculate itself. It controls the finiteness of the vectors and rejects zero or parallel edges: the source would divide by zero in normalization. Almost parallel edges in relative tolerance 1e-12 are rejected as numerically degenerate.

If NB is not written out, the source may round NB to 0 and subsequently divide by zero for a very narrow window. The definition prior to writing checks native calculation for the positive ratio and requires its result >=1; if necessary, it is necessary to explicitly set NB. NB cannot be replaced by an unwritten constant default 50: at other vector lengths, such default would be incorrect.

Fortran array reader allows a shorter field and copies the last folder to the remaining positions; the definition requires explicitly three components so that geometry is not hiddenly supplemented.

## Common defaults and modes

The definition contains CONTROL, TAU, MODE, STRCONST, ENERGY, TASK and SITES. Uses `kkrgen`, `mpi=False`: routines write shared fixed filenames and do not have their own job distribution among the wounds. Numerical or MPI run was not verified here.

ENERGY uses explicit GRID=1, NE=1, ImE=0, EMIN=-0.2 Ry. These are not spatial graph resolution parameters, but initialization options: INIT_MOD_ENERGY and for task with power loop switched off requires NEMAX>0 (`init_mod_energy.f:192–197`). Higher NE does not improve spatial resolution; the definition draws attention to it. EMIN corresponds to native default. Unwritten EMAX=None delegate EF and SITES.NL=None potential context. Listed defaults can be changed; there is no reason to take the zero initial NETAB instead of the applicable initialization.

DMFT and LDAU explicitly prohibits GEN (`gen.f:211–222`); the definition rejects them. ASA or FULLPOT is generally not prohibited: they have separate routines with different outputs. The mode is determined from the potential, not by an added TASK switch. Source does not introduce a general requirement IREL=3.

## Actual behaviour and problems SPR-KKR

1. **FULLPOT ends in experimental PLY routine.** FPPLOT2DPR unconditionally calls PLOT2DPR_FP_PLY (`fpplot2dpr.f:71`). The one with active potential/density task writes PLY files and calls STOP_REGULAR with text `run of test version completed` (`plot2dpr_fp_ply.f:926`). The following regular grid FPPLOT2DPR shall not be carried out during this journey. It is therefore not possible to promise its files `_V`, `_B`, `_RHOCHR` as a normal output of this revision.
2. **Selective flags in experimental FULLPOT path do not distinguish individual groups.** PLOT2DPR(1) and (2) uses PLY routine only for common PLOT_V_RHO (`:101`). Additional potential and density calculations are not conditional on individual flags. PLOT2DPOT/PLOT2DRHO thus do not have the same selective behaviour as the following older branch. Specific PLY files must be evaluated, not the documentation of the older grid writer.
3. **ASA ignores density selection.** PLOT2DPR_ASA does not use PLOT2DPR flags and import VT, not density field. Writes `PLY-V-ASA-*.ply` and `POT-ASA-LINES.ply`. I PLOT2DRHO thus produces potential files in ASA. The definition warns in this variant; because input itself does not recognize the content of potential, the warning describes both branches.
4. **Older FULLPOT branch contains another premature STOP.** After writing `FULLPOT-pot.ply` follows simple STOP (`fpplot2dpr.f:671`); if an end in PLY routine would be deactivated, potential task would still not reach later tables and density. This is the second obstacle to the source, not the path now normally reached.
5. **Uncontrolled degenerate geometry.** Orthogonalisation of parallel edges normalizes zero rest; zero edges and NA/NB=0 lead to further division by zero. The definition protects explicit input by geometric validation.
6. **Strong output names.** PLY routines use names independent of CONTROL.DATASET. Repeated runs in the same directory can conflict; the insulation of these outputs cannot be derived from DATASET. Specialised output parser or change of output names is not part of the patch.

Outputs are based on read radial data and interpolation of geometric window. The number of points and details of the experimental FULLPOT model have other fixed constants that are not issued for user parameters. Numerical accuracy PLY output not verified.

## Verification

16 new cases: discovery, omitted geometric defaults, file/text round trip, pickle/copy, custom window, selective keywords including autodetection when loaded, warning to ASA density, zero/even vectors, orthogonalization, rounding NB to zero, positive NA/NB, vector dimension and prohibited correlation.

Whole set `src/ase2sprkkr/input_parameters/test`: **124 passed**, to default revision plus this separate patch (without previous patches). One expected loading warning PLOT2DRHO reminds the behaviour ASA. `git diff --check` and the reverse application of the patch has passed. SPR-KKR has not been started; tests verify the definition of input.
