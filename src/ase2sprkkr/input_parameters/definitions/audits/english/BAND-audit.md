# 46. BAND in KKRSPEC excl. separate definitions and audits

BAND already exists as TASK option in ARPES, now has a separate definition InputParameters.create("BAND"). Executable kkrspec. Preliminary **works but has major caveats***: complex Bloch roots and their filter cannot be confused with normal BSF or EKREL output. The Numerical Solver hasn't been launched.

## Meaning native calculation

SPEC_RSLAB recognizes BAND and sets IBLOCH=0. SPEC_BLOCH in spec_band.f solves a generalized problem and creates a complex kz z -i*log(eigenvalue)/BRUSEP, with further correction AKCP. It adopts roots according to abs(Im(kz))<EIG_IME_TOL. EIG_IME Therefore, there is no energy in Ry or linewidth in eV; it is a cutoff imaginary part of Wavector in native measuring convention.

SPEC_INPUT at NT=1 prints the label PROJECTED BAND STRUCTURE and at NT>1 NORMAL (along kz). However, SPEC_BLOCH itself stores real kz in K_EIG for NT=1 and KPARA for NT>1. For the interpretation of k-columns, these are actually written to the input variables, not taking the general label from the reader. SPEC_VBRUN subtracts EFERM from E_EIG at the end.

## Separate definition

Shares structure and preparation with ARPES but sets TASK BAND, CONTROL ADSI BAND and adds SPEC_BAND:

- NEIGMAX: positive capacity of stored roots, default 100000 Native ends STOP when exceeded.
- EIG_IME: final positive cutoff, explicit 0.01. Zero would not accept exact real root at strict native inequality.
- CALC_TRANS_PROB: optional attribute of calculation TRLEED transmission properties; default false is not written and corresponds to MOD_SPEC_BAND default.

Default K1={0,0},NK1=1 is written and selected by explicit moment point k_parallel=0. Native has a valid single point to the branch. The path can be set K1 and NK1>1 or K2/NK2 if applicable. These are shift vectors added to KA as in SPEC_INPUT, not endpoints BLOCHSF.

The validation requires this explicit moment procedure, positive numbers, final vectors and their numbers, one energy mesh without SPLITSS and the final non-negative workfunction/mute. Refuses K3/K4/NK3/NK4 which reader does not read. INPVER=0 is governed by the external old input.inp, the content of which does not Verifyme here.

## Defaults: explicit deviations and omitted behaviour

In BAND definitions are **written** IMV_INI_EV=0 and IMV_FIN_EV=0 eV. For non-absorption band procedures are more suitable than taken photoemission 0.05/2 eV, but for complex/disordered band workflow the user can choose non-zero control. Nejde o elektronovou teplotu.

EIG_IME=0.01 is **written**. Native starts at 0.01 but at single EIG_IME and VIL>1e−4 sets 9999, thus practically shuts down the filter. Python therefore does not promise unwritten default 0.01 and does not leave the filter dependent hidden when changing the control.

NEIGMAX=100000 is written to the input and corresponds to native initialisation. Zero CALC_TRANS_PROB is omitted from the input and retained native false. CRYS_VECS=False is omitted from the input and retains the use of primitive crystal vectors; presence of correctly named key native turns off primitive vectors. Original ARPES CRYS_VEC is not read and corrected in the new definition.

Listed NE=300, EMINEV=−8 and EMAXEV=+5 eV against EF, EWORK_EV=4.2 eV, NL=4, N_LAYER=50, NLAT_G_VEC=57 and N_LAYDBL={10,10} remained practical choices. For belts comprising both occupied and vacant states, +5 eV has a sense as opposed to valence AIPES. The number of points and cutoffs must be converged; native smaller initialization numbers are not automatically better.

SPEC_PH is left as part of shared native spectroscopic preparation. It cannot be assumed that BAND is the photoemission intensity or that EPHOT is the shift of all energy bands: SPEC_VBRUN for IBLOCH=0 sets EHIGH=ELOW before FINALSTATE. This audit does not certify any dependence of auxiliary preparation on photon/angle elections.

## Specific native errors and ALE

1. spec_band.f:70–80 at an inappropriate ratio EIGA/EIGB writes SPEC_BLOCH: not converted but does not perform CYCLE. EIGENVAL remains zero and follows CDLOG(EIG). The failure thus reported cannot be ignored; it can create an invalid complex kz. Python definition to does not fix.
2. The imaginary kz filter and non-absorption default is not proof of the physical correctness or existence of all belts. Numerical tolerance, energy mesh, surface direction and reciprocal basis affect the roots taken.
3. CALC_TRANS_PROB is native available but transmission amplitudes and their normalization have not been numerically verified. Default false preserved.
4. Native BAND in SPEC_VBRUN enforces ATA=1; disorder/CPA result therefore requires own interpretation and is not generally certified by this definition.
5. mpi=False is a conservative setting. Native has its own parallel collection of band data over rank files; this string has not been numerically tested. There is no evidence of error of all MPI calculations.
6. Surface structure, workflow, FULLPOT, many-body and thermal branches are not completely certified by this input-only audit.

## Tests

10 new tests; whole set InputParameters **118 passed**, two existing NumPy DeprecationWarning messages. Discovery verification, band defaults, zero control, filter, capacity and unwritten transmission option, round-trip/pickle, executive, path moment and transmission flag; invalid cutoffs, numbers, ignored vectors, negative control and SPLITSS.

Command: `PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 ../test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q`.

The sourcee: spec_rslab.f:160–225; spec_inputs.f:1397–1418,1420–1463,2223–2237; spec_band.f:66–100; spec_vbrun.f:203,224–230,566,597,786–810; mod_spec_commons.f:232–242. ASE baseline 267512c4b8ccbbfd37af7ede898febc066b7e3ee. Single patch only for BAND; without dependence on previous patches.
