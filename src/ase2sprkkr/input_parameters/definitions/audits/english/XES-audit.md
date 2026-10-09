# XES new definition and audits

Delivered SPR-KKR 9.7 compared to fetched branch ASE2SPRKKR refactor, commit 267512c4b8ccbbfd37af7ede898febc066b7e3ee. XES does not have its own definition. Patch 0024-xes.patch adds xes.py and test_xes.py separately against this base.

## Interface

INITVAR recognizes TASK.XES and sets X-ray workflow IXRAY=1, ICORE=0, RDTAUDIR=true and INITELOOP=true. GEN:1378–1385 calls XRAYSPEC which internally renames XES to RXES. CHECK_INPUT:43–47 requires IREL=3. The definition therefore rejects other MODE.

XRAYSPEC:376–380 requires ENERGY.NE and its equality of points from the first TAU contour. The definition itself does not read the content TAU directory: its mesh, number of points, types and potential must correspond to the actual calculation. When available TAU is read, otherwise the main workflow can prepare the loop; the numerical link has not been started here.

TASK exhibits IT, CL, ME, OUTPUT, POLOBS, DIPOL, GFCONVRH, WRREGIRR, NONDIPE2, NONDIPM1, TAUCORE and Euler angles FRAMETET/FRAMEPHI/FRAMEGAM. This is not an interface ARPES or VBPES. Frame angles rotate a common coordinate frame, not a simple single vector of an electric field; init_mrot_frame.f:48 form a matrix GETMROT(-gamma,-theta,-phi).

## Defaults and writing

| Option | Native effect |Definitions|
| --- | --- | --- |
| IT | 1 | 1, omitted |
| CL |2P from common reader| 2P, omitted |
| ME | ADA | ADA, omitted |
| OUTPUT |BARN; normal output in Mb| BARN, omitted |
| POLOBS | S | S, omitted |
| DIPOL |0; for NAB it is overwritten to 1| 0, omitted |
| GFCONVRH |false in absence of keyword| false, omitted |
| WRREGIRR |DATA lists true but SECTION_FIND_KEYWORD overwrites them by presence| false, omitted |
| TAUCORE |Two zeros, hodnoty<1e-6 take over table WCOREHOLE| None, omitted |
|Euler angles| 0/0/0 |None, omitted; absence does not switch CHANGE_FRAME|
| GRID/NE |XES does not have its own initialization default; NE=0 leads to invalid capacity|explicitly3/300|
| ImE | 0,01 Ry |explicitly0,01 Ry|
| EMIN | -0,2 Ry |explicitly -0,2 Ry|
| EMAX |XES explicitly overwrites to EF| None, omitted |

GRID3/300 is an explicit choice of the energy scan of occupied states, not an unwritten statement pretending to conform to native default. Convergence EMIN/NE/ImE not verified. Native NE=0 is not more appropriate; other baseline values are preserved. Deleted EMAX, effective WRREGIRR=false and Euler angles correspond to actual native behaviour.

The definition limits n<=5 because SHELL contains only K,L,M,N,O, and l<n. IT must be positive; its upper limit and the existence of core shells require potential. NE>=2 protects GRID3 from uninitialized ESTEP at one point. Split and two paths are rejected due to the use of the first TAU contour. ME=NAB s DIPOL=0 fires overwrite alerts; NONDIPE2 takes precedence over NONDIPM1.

## Major caveats and native errors

1. **MPI spectrum does not start at all.**GEN:1355–1368 in the list of continuing tasks does not have XES; after the energy loop STOP_REGULAR ends running before XRAYSPEC. The definition has mpi=False. A similar risk relates to an existing XAS definition with mpi=True; a complete XAS audit will be a separate step.
2. **RXES shares its own rate calculation with RXAS.** The only occurrence RXES outside the comment is when renaming task. The calculation and output branch distinguishes XMO/XRS against all others, so XES enters into the common absorption-rate branch. Therefore, the result cannot be declared as a verified emission spectrum without physical comparison with the reference. It does not have a simple proxy STOP as FORCETEN so that task has not been skipped but will have a significant ALE in a joint audit.
3. **The output is still marked XAS.** The .rat name and text SPEC use the same core shell convention as XAS; the text states XAS spectrum. XAS and XES for the same DATASET, type and shell may write the same file.
4. **MESPLIT is not a switch:** after loading the keyword, K_ME_SPLIT always sets true (297–298). The definition does not expose it as an effective choice. DYNAMICAL is read only for RXAS/XMO, therefore it is not offered for XES.
5. **TAUCORE is read but its use for RXES rate is not documented.** Values are transferred from eV to Ry or taken WCOREHOLE. The subsequent use is in extended XMO/XRS branches; the standard matrix-element call sets CECOR=0. Parameter is not designated as effective emission Broadening management. Change TAUCORE does not have to change the spectrum itself RXES.
6. **ME is not proof of three full Gauge calculations.** MEFORM is loaded; NAB forces DIPOL=1, GRV prepares radial derivations. However, the actual normal matrix elements is calculated ME_ALF_ALF, without delivery MEFORM. Only acceptance ADA/NAB/GRV can not be confused with verified independent implementation of all three forms.
7. **Non-dipole variants depend on formalism.** ME_ALF_ALF_NONDIP is called only for .NOT.FULLPOT and DIPOL=0 (805–812). For FULLPOT the same flags do not automatically guarantee the same contributions. Potential formalism is not read pure InputParameters validator.
8. **Only one atomic type:** NTXRAY=1 and IT is loaded as a single number; it is not a list of types. The output has three polarization records for each core state and energy, not a single intensity column. POLOBS=S gives (+,-,z), C gives (x,y,z).
9. **Standing XAS has an ineffective default EMAX in Python:** ENERGY(defaults={"EMAX":4.0,...}) does not use this key. XAS natively sets EMAX=4 Ry, but the declaration does not issue the limit as a true member. The new XES module does not take over this misleading design. Another extension XAS is separate.

## Tests

16 new cases: serialization, reload, pickle; omitted defaults; type/core shell/polarisation/width/frame; relativity limitation, mesh and split contour; invalid IT/DIPOL; physical shells and limit n=5; NAB override. Full input_parameter set: **124 passed**. Numerical SPR-KKR has not been started, emissions or .rat data have not been verified. Patch passed git diff --check and re-application control.
