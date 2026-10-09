# TUTORIAL Total audit SPR-KKR 9.7 / ASE2SPRKKR

This document replaces sub- TUTORIAL audits. The whole definition has 14 options TUT-TASK. Implementation prepares and validates input; it is not a repair of Fortran solver or a confirmation of the numerical correctness of all branches.

## Checking the previous file
Before the last extension tutorial.py actually included POT-XC and CORE-ALL, including info and validation FULLPOT. For a clear download, it is now included tutorial-complete.py; its content is identical to the current definitions/tutorial.py and is installed under the name tutorial.py. The previous implementation division after branches is completed.

## All Options
| TUT-TASK |Calculation|Status / main conditions|
|---|---|---|
| YLM |Real spherical harmonic as spatial area; positive/negative part separately|Fully meaningful for instruction; fixed output names|
| JLX |Spherical j_l(x), n_l(x)| Fully meaningful; LMAX<=4, NX>=2 |
| JLX-ASYMPT |The same against small/large x asymptote| Fully meaningful; SMALL vynucuje XMIN=.001 |
| CHR-POT |Saved V/B, charge/spin density, FULLPOT shape functions|Fully meaningful for diagnosis of available data|
| POT-DECOMP |Nuclear, Hartree and MJW XC/X potential|Major caveats: ASA, original XC VBH/MJW/VWN, NONMAG bug|
| POT-XC |VBH/MJW/VWN/PBE at the same density|Major caveats: ASA, NONMAG bug, fixed graph range|
| RWF-L |Regular radial functions across l at EMIN|Major caveats; additional matching graph is broken|
| RWF-E |Regular radial functions across energy at specified l|Major caveats: scalar relativity, ASA, fair part and normalization|
| CORE-ALL |Non-magnetic core states of type 1|Major caveats; without core shells broken branch|
| DOS-FEG | V=.0001 Ry, B=0, single-site DOS |Meaningful demonstration; STANDARD KKR, not full lump free-electric integration|
| DOS-SS |Single-site DOS original potential|Meaningful in STANDARD KKR; CPA not certified|
| DOS-CLU | DOS konfigurace clusteru |Meaningful with explicit TAU cluster; name does not switch on|
| DOS-BULK |DOS periodic system|Meaningful of periodic BZ integration|
| GEN |General GEN branch a DOS/observables|Meaningful, result determines configuration TAU/potential|

None of these 14 options are generally off placeholder. The breakdown relates to this task; the joint audit of all tasks is not yet created.

## Native Dispatch and inputs
initvar.f:700–718 requires TUT-TASK; for DOS* and GEN turns on the energy loop. check_input.f:167–172 requires KKRGEN. tutorial.f:53–184 contains the entire dispatch table, upon completion STOP_REGULAR (normal termination). The exact canonical values are declared in Python; native prefix variants are not added as another tasks.

YLM reads L/M (2/0), silently cuts M. JLX reads LMAX/NX/XMIN/XMAX (2/300/.001/5), cuts LMAX to4, fixes the limits. JLX-ASYMPT plus X-RANGE SMALL/LARGE. The validation rejects invalid limits/index instead of silent corrections.

RWF reads TASK.ITSEL (1), NORMALISATION (1), ENERGY.EMIN/EMAX/NE (1e-7 Ry/1 Ry/10). NORMALISATION=1 integrally normalizes ** the displayed real part** radial functions; any other native string preserves MST amplitude, the definition uses canonical MST. RWF-E plus TASK.L (2, native clamp to NL-1). Existence of type ITSEL and actual NL require potential; L we control against SITES.NL if it is explicit.

CHR-POT sets the whole PLOTPRS=true and calls FPPLOTPRS. Routine fpplotprs.f branch outputs according to FULLPOT and available fields; it does not have other own reading input parameters. This selection does not need ASA restrictions or restrictions on original XC.

DOS-FEG overwrites VT=.0001 and BT=0, DOS-SS maintains VT/BT; both sets IBZINT=1, NPTMIN=NPTMAX=0. V tau_drive.f:88–102 is tau=t used only below STANDARD-KKR; therefore, the definition of this mode explicitly writes and rejects the cluster for these two options. It cannot be confused with a normal lump FEG integration. DOS-CLU/BULK/GEN only call GEN s DOSCORSYS=LOC; do not switch geometry themselves. DOS-CLU requires TAU.CLUSTER/MOL/BZINT=CLUSTER and CLURAD or NSHLCLU; DOS-BULK rejects cluster. GEN maintains the set geometry. GEN for TASK=TUTORIAL turns off core calculation (gen.f:394) and writes CALCDOS_WRITE.

## Defaults: written versus absent
Default option YLM is explicitly written Python option; native without TUT-TASK stops. Mathematical defaults agrees with native.

ENERGY defaults is changed according to the current TUT-TASK unless the user overwrites them:
- mathematics/plots/core: GRID0, NE1, ImE0; limits -.2..1 Ry are only startup values, these branches do not use them;
- RWF: GRID3, NE10, EMIN1e-7, EMAX1, ImE0; Routine itself shall reconstitute EPATH(3), with zero imaginary part and ignore input GRID/ImE/eV limits. Definition rejects relative eV limits, more contours/SPLITSS and NE1;
- DOS/GEN: GRID3, NE300, EMIN-.2, EMAX1, ImE.01 Ry. Everything is explicitly writes; it is a practical spectral option, not native default contour. EMAX1 allows a vacant part; the alternative is a window up to EF according to the purpose of the calculation.

MODE explicitly SREL for RWF, FREL for others (native FREL is absence NREL/SREL/SP-SREL; token FREL does not change this result). TAU.KKRMODE=STANDARD explicitly in all default variants; the user may change it except DOS-FEG/SS. Cluster v native init_mod_kmesh sets REAL-SPACE according to IBZINT0. BZINT defaults POINTS, s CLUSTER/MOL CLUSTER; this avoids the original wrapper problem callable None transferred to the string NONE.

SCF section contains only FULLPOT flag, at default is absent. This is not a guarantee ASA: native takes over potential format. There is no false SCF.VXC override: this value for POT-DECOMP reads XC-POT from potential (fallback VWN v potio.f:963), POT-XC internally switches it.

## RWF native errors and limits
TUTORIAL_RWF (tutorial.f:1422–1760) uses WAVFUN_READ_SRA, therefore the definition allows only NREL/SREL/SP-SREL, default SREL; FREL rejects. IRTOP is set JRWS(IM) while FULLPOT reader expects JRCRI(IM): FULLPOT is excluded from the audited range. Numerical energy close to zero changes internally to 1e-5 Ry.

RWF-L the main graph creates l=0..NL-1 at the first energy. Another matching graph has fundamental errors: YMATCH(I,1)=YY(I,0,0), with YY(:,0,:) remaining zero; L following the loop corresponds to the last channel. For external radius, YMATCH(:,2) but YMATCH(:,1) is not to be completed. The matching graph therefore contains zero internal wavefunction and uninitialized external data. The definition itself does not correct this native problem.

RWF normalization has 1/sqrt(norm) without guard for zero standard. The output draws DREAL(ZG), not complex wavefunction or the whole relativistic spinor amplitude. INDEX0 is a zero reference curve in the main graph. The routine translates ETAB/WETAB, uses NEMAX from initialization and requires consistent one mesh; the definition holds the same NE at startup and reread. The validity of potential data/ITSEL cannot be guaranteed without their loading.

## Parallelism and extent of verification
mpi=False for the whole task. GEN MPI guard does not allow TUTORIAL, mathematical and fence branches use common files. No numerical SPR-KKR calculation was performed.

40 new tutorial regression cases; full input_parameters set: 148 passed, two original NumPy warnings. Covered by all14 canonical options roundtrip, pickle RWF, change dependent defaults, geometry DOS, mesh, explicit FULLPOT, relativity, energy references and mathematical indices. Tests verify Python definition, not accuracy native outputs.

## Detailed findings of potential/core branches
## What really counts
For each type ITBOT..ITTOP and radius R up to JRWS creates five curves: V_tot taken from VT, V_nuc=-2*Z/r, residual V_H=V_tot-V_nuc-V_xc(previous), MJW V_xc and MJW exchange-only V_x. The original V_xc shall be re-evaluated from the stored RHOCHR/RHOSPN function marked in the potential. For magnetic types creates analog B curves from the spin difference. Output is Grace/XMGR file created XMGRHEAD with DATSET0, key Vdecomp and TXT_T, radius in a.u. and potential in Ry. It does not perform SCF or DOS.

Attention: Drawn MJW V_xc is not always V_xc(previous) deducted during the construction of the Hartree residue. For input VBH/VWN therefore V_nuc+V_H+drawn V_xc does not generally reproduce V_tot. V_x is a separate exchange contribution, not another member to be added to V_xc.

## Actual restrictions
- FULLPOT: tutorial.f:633 performs RETURN without output. The new expert of the duty section SCF.FULLPOT allows the explicit request to be rejected. FULLPOT contained only in POTFIL without reading the potential of the Python validator.
- Original SCFVXC must be VBH, MJW or VWN (r.671–693). PBE leads to STOP, other unsupported functionals also. The value is read from XC-POT in potio.f:963; the fallback there is VWN. SCF.VXC is read in scf.f:384, this tutorial branch SCF does not start. Therefore, it is not exposed to false override or invented default functionality.
- The Hartree item is a residue against the saved total potential, a nuclear member and a reconstructed XC, not an independent Poisson equation. It requires consistent density/potential data.
- NONMAG is a shared status field: at abs(B) below 1e-4 it is set true (r.722–723), but in the following type it is not returned to false. A combined system where nonmagnetic type precedes magnetic may lose B graphs. Already the initial NONMAG=true B will suppress the graphs. It's a native error, irreparable by itself input definition.
- Serial settings TUTORIAL preserved; common output files are not certified for MPI.

Category for future joint audit: works but has a large ALE. FULLPOT path does not provide a result; it is not a generally disabled placeholder.



## What compares
The hard encoded list VXCTAB (r.847) contains VBH, MJW, VWN (LDA) and PBE (GGA). For each type ITBOT..ITTOP evaluates all four functionals on the same saved RHOCHR/RHOSPN. For PBE prepares spin density, radial derivation, Laplacian and grad rho . grad rho' (lines 934–987). Output Vxc=(Vup+Vdown)/2 and Bxc=(Vup-Vdown)/2 writes to Grace/XMGR to files with DATSET0, key Vxc and TXT_T; radius in a.u., potentials in Ry.

It's not about four separate self-consistent solutions. Default XC-POT potential is not limited to VBH/MJW/VWN as for POT-DECOMP: here SCFVXC internally switches between fixed four functions. We do not introduce input SCF.VXC because this branch does not manage the list. The quality of the result depends on the stored density and its radial derivations.

## Restrictions / Mistakes native Implementation
- FULLPOT leads to a silent RETURN (r.867). The verifier rejects explicit SCF.FULLPOT; the full-potential format itself POTFIL shall not be verified without reading potential.
- NONMAG is set at a small Bxc true (lines 1028–1029), but for the next type is not returned false. The non-magnetic type before magnetic can therefore suppress Bxc charts of the later type; the problem is similar POT-DECOMP.
- The Vxc graph range is converted to -15 after the calculation..0 Ry (r.1031–1032). Data outside the range are written but can be hidden in the normal display. It's not a numerical cutoff.
- SCFVXC remains after the calculation PBE, the original value is not renewed. TUTORIAL then normally ends STOP_REGULAR, therefore in standard standalone run does not follow SCF with unexpected function. In future insertion of the routine into longer workflow it would be a significant side effect.
- PBE preparation highlights positive radial density derivations but continues; expressions GDGAG are compiled as -n'*n'. This is not a confirmation of correctness for any density profile.
- MPI remains off, parallel outputs are not certified.

Category for future joint audit: works but has a large ALE. The placeholder isn't off.



## The Real Meaning
Despite the name ALL, it uses IT=1 (lines 1199), not all atomic types. BT is reset for all types and NONMAG is set true (r.1177–1178). Recalculates nonmagnetic core states in spherical potential VT of the first type, not original spin split. Z single shells selects m_j=1/2; saves the first radial component for each selected state from the calculation CORE.

Selection shells: table n,l and estimate according to Z; if NCORT(1) !=0, count shells shall be derived from cumulative filled shells corresponding to NCORT. If no match exists, native STOP. It is not freely adjustable CL selector.

Grace/XMGR output contains VT(r), core level energy and radial g/f component scaled and shifted to level energy: E+0.2*g(r), E+0.2*f(r). Therefore the vertical position of the curve is not the original amplitude wavefunction. Radius is in Bohr, energy graph in Ry; console also lists eV, standards and other quantum numbers. Energy is not deducted from EF here. The file consists of XMGRHEAD from DATSET0, core-all keys and TXT_T(1).

## Native errors / restrictions
- FULLPOT leads to STOP_MESSAGE (lines 1153–1156) rather than to a supported computing branch. The verifier rejects explicit SCF.FULLPOT; the format itself POTFIL shall not be verified without reading potential.
- Without core shells, e.g. Z<=2 and NCORT=0, not NNKCOR initialized: it is set up in the found core state cycle, but subsequently used for the number of curves and other loops. At the end of the empty shell cycle, an additional MJ/LCXRAY cycle is used. This branch is broken. Definition without reading potential does not guarantee empty core states.
- When writing VT, the length JRWS(IT)-IVBOT+1 is used instead JRWS(IM)-IVBOT+1 (r.1366–1367), while the data and radius are indexed IM=IMT(1). For the first type with a mesh index other than 1, the length may be incorrect; reading outside the valid radial area is at risk at different limits.
- BT, NONMAG, NCXRAY and LCXRAY are mutated. V normal standalone tutorial run routine ends STOP_REGULAR; these mutations are not saved as new potential.
- MPI remains off; serial outputs are not certified for parallel use.

Category for future joint audit: works but has a large ALE for ASA type potential 1 with blank consistent core shells; without core shells broken branch. It's not a generally off placeholder.


## Delivery
0053-tutorial-complete.patch is the total cumulative baseline patch 267512c4, replaced by 0049–0052. Do not use old tutorial patche simultaneously. Replace tutorial.py and test_tutorial.py with the entire files after an old patch has already been applied.
