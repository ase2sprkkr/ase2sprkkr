# Consolidated audit of SPR-KKR 9.7 and ASE2SPRKKR

Status as of 4 October 2026. Sources: the supplied **PUB9.7_260420_18h13** code and the preceding individual audits of InputParameters definitions. This report consolidates the findings according to the actual reachability of each calculation and the meaning of its output. It is not a numerical certification of SPR-KKR.

## Scope and meaning of the categories

The audit covers the tasks, aliases, and selectors reviewed so far in INITVAR/GEN/RSP/KKRSPEC, including the complete TUTORIAL task. It is not a comprehensive audit of every internal algorithm, correlation method, or experimental combination of switches.

| Category |Meaning|
|---|---|
| **Fully meaningful** | The reachable implementation matches the exposed purpose within the stated scope; no fundamental obstacle to that workflow was found. Ordinary formalism restrictions or convergence requirements do not by themselves constitute major caveats. |
| **Caveats — 2A: mostly functional** | The main audited path is meaningful; it requires a suitable state or input, or restriction to a working branch. |
| **Caveats — 2B: mostly nonfunctional for the stated purpose** | The main declared or explicitly identified branch is unfinished, produces defective data, or has an unsubstantiated physical interpretation. Supplying the correct state alone is insufficient. |
| **Broken** | A demonstrated defect in reachability, initialization, indexing, or data flow blocks the normal purpose or corrupts the result. An input definition does not fix it. |
| **Disabled** | The relevant path is explicitly stopped or replaced with an unavailable routine. The combined list appears at the very end. |

The classification applies to the supplied SPR-KKR revision. For a task with distinct branches, the scope is stated explicitly. Valid input, a passing pytest run, or normal program termination does not establish a correct physical result.

## 1. These tasks are fully meaningful

|Task / Branch|Purpose and audited scope|Important support boundaries|
|---|---|---|
| **SCF** |Self-consistency calculation; corrected options BROYDEN2, ANDERSON, TCHEBY.|PULAY is not valid option of this native reader. TCHEBY is not for FULLPOT. The export and special correlation workflow has not been generally certified.|
| **DOS** |Standard DOS over generic GEN loop; native internal TASK is NONE.|GRID, numbers and windows shall correspond to the required spectrum. Advanced DOSREP branches have specific MPI risks and are not automatically covered.|
| **BSF** |Bloch spectral function; energy scan EK and single energy cut KK.|KK is not automatically to EF. Existing EK default EMAX=-1 Ry is a suspect explicit interval; see defaults below. All path mapping and 2D transformations have not been fully verified.|
| **EKREL** |E(k) and Fermi surface for arranged systems.|CPA is rejected natively. TASK.EMIN/EMAX/NE controls its own search, not ENERGY. The preset to the path depends on the grid; ESLICE=EF may be outside the search interval.|
| **WFPLOT** |Diagnostic charts of wave functions BAND/CORE in ASA, FREL.|It uses the first energy point; availability shells and atomic type must determine potential. FULLPOT and those correlations are not supported.|
| **CLXPS** |Core-level XPS in ASA/FREL with implemented ADA branches.|Only implemented forms ME; the final state must be above EF. Geometry has actual sufixes 0 or continuous 1..10. Capacity ENERGY.NE must be sufficient to core states.|
| **CHI** |Audited static response in ASA/FREL, POINTS, serial.|NCHIITER is natively always 1; it is not an iterative solve controlled by this option. Local CHIORBPOL is not identical to general MODE.OP.|
| **POSANITAU** |Positron annihilation integration through own TAU path in ASA/FREL.|The suitability of the energy limit EF and the zero partial speed shall be assessed; 1/RATET is not protected from zero. This is not a substitute for preparation POSANI.|
| **JORB** |Orbital current density, FREL, serial.|The potential must provide the necessary angular representation. Auxiliary fort.* are not a clear main output; MPI completion does not work.|
| **PLOTPOT** |Selective chart of saved potential, in ASA and FULLPOT.|Does not perform SCF. The data and their mesh must be present in potential.|
| **PLOTSFN** |Graph of stored radial coefficients shape functions in FULLPOT.|Does not create SFN or 3D cell shape; FULLPOT flag does not convert ASA potential. In ASA, the SFN plot is not produced.|
| **TUTORIAL: YLM, JLX, JLX-ASYMPT** |Mathematical examples of harmonics, Bessel functions and asymptotics.|L/M and limits are validated; even mathematical branches undergo common initialization and need usable input environment.|
| **TUTORIAL: DOS-FEG, DOS-SS, DOS-CLU, DOS-BULK, GEN** |Samples DOS / normal GEN workflow in the corresponding model.|FEG and SS are single-site STANDARD. CLU requires actual cluster geometry and CLURAD/NSHL; the name of tutorial alone does not create it. BULK requires periodic integration.|

For BSF needs to be corrected or intentionally confirmed default window before use without its own parameters. U POSANITAU includes an unprotected zero partial contribution among specific defective cases; does not mean the fact that each non-zero total annihilation integration is broken.

## 2. Tasks with caveats — classification by practical usability

### 2A. Mostly functional: require the correct state, input, or restriction to a working branch

V for this group the main audited calculation has a meaningful reachable path. Typically, there is a need for suitable potential, compatible TAU/DOS, correct magnetic geometry or use of only implemented formalism. It's not the same as the unfinished main solver. Some tasks also have defective secondary output or alternative branch; these are explicitly set aside to 2B.

Classification to 2A is not confirmation of the numerical correctness of all combinations. ©The correct status here means the conditions of a specific physical workflow, not any tuning of parameters that only hides a programming error.

| Task / rodina |Condition of use or restriction|Practical Conclusion|
|---|---|---|
| **MEPLOT** |A photon scan with more electron energy rewrites curves. The presence NEPHOT=1 itself switches the reader to the range branch; SEFIN receives the initial ERYD.|New definition does not write single-photon default NEPHOT and forbids current double scan. Meaning SEFIN for the final state still requires physical verification.|
| **MECHECK** |NMEFORM is overwritten, the next test is complete-time disabled and the summary loop has an incorrect upper limit.|Sub-total ADA/NAB comparisons make sense; all checks described cannot be promised GRV.|
| **SOCPAR** |Test ABS(... )<0 never reverses the mole. EXCME is missing Ry→eV conversion, although the axis says eV.|SOC part and exchange part must be interpreted separately. Exchange output has documented non-compliance of units; correct Fortran, not input default.|
| **PSHIFT** |Phase shifts exist, but resonance fracturing may be used as an unfound resonance; interpolation does not control actual pi/2 exceedance.|Phases are not automatically broken; derived resonance numbers are not reliable without repair. NR non-uniform cutoffs and high relativistic legends have other risks.|
| **PLOTPRS** |Common task bypasses the error of the separate PLOTRHO, but the spin and orbital curve points to the same dataset and the panel limits ignore the orbital density.|Potential, charge and SFN outputs make sense; resolution of magnetic density in the graph requires correction of Writer.|
| **STONER** |Kernel is fixed VWN; small DOS is replaced by nonphysical unit contributions; some normalization is separated without protection.|The result is not automatically the kernel used SCF of the functional. It needs suitable cutoff and nonzero standardization variables.|
| **VBPES, VBUPS, VBXPS** |GEN implementation is common; VBUPS does not switch to UV photon. DEFEP0/EWORK are only diagnostic and custom geometry can be overwritten. TAU and caches are not verified.|UV specify the energy. This family is not identical to KKRSPEC ARPES/AIPES.|
| **APS** |Reads old DOS format without capacity control. Broadening expands more fields and depends on absolute energy.|Without Broadening, there is a compatible path DOS. GAMMA is not normal constant Lorentz width; positive GAMMA requires correction or specific verification.|
| **AES** |Only NQ=1; cache ME/MEi does not control origin. NE of TAU may exceed the allocation. F shells may exceed the table of names.|The definition excludes f shells and NEME<3; compatible uniform growing TAU mesh and pure cache remain necessary.|
| **FMAG** |Spin/orbital form normalization divides by moment at q=0 without protection. Network q is fixed.|Zero moments can give NN/Inf even in a otherwise correct system. The normalized curves cannot be recounted as the original moments.|
| **JALF** |The main current file has a purpose, but the auxiliary diagnostic mixes types and meshes; at NL<3 exceeds the field.|The definition protects explicit cutoff. Vector potential does not count as standard, although the code contains an unachieved branch.|
| **TORQUE** |Needs recorded moments and appropriate torque axis; the normal unrecorded potential may end up STOP. Mixed cast uses the first type.|The validation of input without potential cannot guarantee an appropriate magnetic geometry or meaning WGT_OCC.|
| **XAS** |The static branch has a purpose but accepted ME names do not show three independent Gauges; non-dipole formalism depends on ASA/FULLPOT. The limits can natively overwrite EF.|Fixed MPI metadata, EMAX and incorrect MECHECK. TAU conformity and effective scan start should be checked outside serialization.|
| **PHAGEN** |Only ASA/VWN, NT=1, without disorder. It exports the original total potential together with density/Coulomb potential after the next step. Spin information is not exported and the radial index may be incorrect.|Export for simple NT1/mesh1 system makes sense. It's not a general converter, it's not even its own phase-shift solver.|
| **JXC** |FREL symmetrical rotation is eliminated by condition OR 1.EQ.1. Cluster SP-SREL path normally does not have prepared scratch files.|Simple periodic collinear workflow makes sense; tensor at equivalent locations and cluster are generally not reliable.|
| **DMI** |Dz is calculated differently from Dx/Dy. The limitation of representatives and couples does not automatically complement the multipleness of the complete cell.|Three columns are not automatically equivalent components of one procedure. Micromagnetic totals require inspection of frames, multiples and radius.|
| **COMPTON** |Native SPLITSS is read from two sections. EMD2D with nonzero PNMAX scales the longitudinal vector twice; even shifted meshes loses the extreme points.|Definition writes both SPLITSS and excludes defective EMD2D/shift combinations. The basic compatible workflow is not marked as broken by this.|
| **POSANI** |Depends on compatible data POSANIPREP; file TAUQQ has NQ/NT error. Core contribution is limited and core lifetime parameterization has a specific fixed shape.|It's not a single universal positron during. Before using prepared data, check their structure and origin.|
| **POSANIPREP** |TAUQQwriter/reader iterates NT instead of NQ: can skip site or cross the field. Old rho cache can skip export waves.|Meaningful preparation in a compatible simple case; general NQ≠NT workflow is broken.|
| **GILBERT** |GEFF divides the spin/total moment from protection. Heat mesh and standardisation need specific physical data.|Basic Damping has a sense of suitable magnetic system; compensated or zero moments are defective cases.|
| **CHILANDAU** |Incorrect finishing condition prints continuous integrals; RDTAUFIL does not prepare the necessary energy content. MPI has allocation/reduction errors.|Serial calculation to assess by the last complete block. Reuse TAU and MPI are not functional general routes.|
| **SIGMA** |Fermi-surface and sea are different branches; the sea is not automatically included. ImE and fixed names sigmaXXXXX.dat affect use.|Do not present the surface itself default as complete all posts. The scope of the definition is audited ASA/FREL/POINTS path.|
| **OPTABS** |NOMEGA is the allocated capacity but the loop uses NEB. TEMPEL<=1 leaves KBT unprepared. NEA is overwritten and non-zero OMEGA1 violates audited pairing.|The definition requires TEMPEL>1, NOMEGA>=NEB and limits the audited window. Task is OPTICS-ABS, not generic OPTICS.|
| **CHIXAS** |It needs external NEF/CHISPN/CHIVV data with the right units, continuous numbering and potential dimensions.|No running with universal physical defaults. Test fixation does not replace measured/calculated input data.|
| **ARPES** |Old TASK aliases, incorrect CRYS_VEC and native IMV_INI in Ry writing another Damping variable.|Use custom AIPES/BAND/SPLEED definition. Warning when validating/ storing an old alias is DeprecationWarning and Python can filter it.|
| **SPLEED** |Shares KKRSPEC but has a different branch and an absolute incident EV. Exhibited is conservative solid geometric workflow.|Not ARPES with a different label. The Damping/unit and surface potential structure shall be appropriate.|
| **AIPES** |IBLOCH=2 sets local/incoherent calculation of all atoms; it is not an ordinary numerical angular integration ARPES.|The difference against ARPES is in the model matrix elements and interference, not just at the output angle.|
| **BAND** |Complex strips use -i LOG(EIG), cutoff is in Im(kz). When failed EIGA/B there is no interruption and may come LOG(0); column labels vary according to geometry.|There is no substitution EKREL for ordinary E(k) Appropriate eigenvalue cutoff is required, column interpretation and solver failure repair.|
| **CLPES** |Core window rewrites the assigned bounds. CPA and some polarization stops running; ECOR in NE=1 path may come from the wrong last type.|Serial fixed geometry and suitable shell have meaning; explicit photon 1000 eV is not a guarantee of crossing the threshold.|
| **TUTORIAL: CHR-POT** |It uses common PLOTPRS, including magnetic graph problems.|It's a graph of read data, not a new SCF calculation.|
| **TUTORIAL: POT-DECOMP, POT-XC** |ASA; NONMAG remains sticky between types. POT-DECOMP breaks down the state with the original Vxc, the shown new Vxc does not have to correspond to the original.|POT-XC compares functionals at the same density, not four SCF solutions. PBE is not supported original XC in the branch POT-DECOMP.|
| **TUTORIAL: RWF-L, RWF-E** |It uses scalar reader, not generic FREL. Normalization has no protection from zero; matching fence RWF-L reads blank data.|Radial data have a diagnostic sense in scalar/ASA range, matching curve RWF-L is broken.|
| **TUTORIAL: CORE-ALL** |Only the first type, zero exchange field and selected output components. Without core states, variables uninitialized may remain; in another mesh index, the error is JRWS.|ALL does not mean all types and full spinor. Suitable first type/mesh/core conditions are not guaranteed by input.|

**Examples of difference:** TORQUE with unrecorded potential is not proof of unfinished torque solver; CHIXAS without proper external susceptibility does not have necessary physical inputs. By contrast, faulty tensor rotations JXC or uninitialized spectral fields are not a matter of appropriate baseline.

### 2B. Mostly nonfunctional for the stated purpose: defects, an unfinished main branch, or unreliable results

I do not consider this group to be normally applicable only after the correct entry has been completed. It is necessary to correct the implementation or verification that the calculation corresponds to the declared purpose at all. U mixed tasks the conclusion refers only to the specified branch; their functional basis remains in 2A.

For XES and TDOBS, the issue is whether a separate implementation of the declared functionality is established. Source inspection reveals this concern, but without a numerical reference comparison it does not prove that every output is wrong. For the other programming defects listed below, a specific cause is demonstrated.

| Task / specific branch | Why correct input is insufficient | Practical conclusion |
|---|---|---|
| **PLOT2DPR / PLOT2DPOT / PLOT2DRHO: FULLPOT** |Unconditional experimental PLY branch writes both groups and end run before older regular gridwriter.|The result is not the promised selective regular tables. NA/NB does not control all PLY sampling.|
| **PLOT2DPR / PLOT2DPOT: ASA** |ASA Writer draws potential; whole PLOT2DPR does not deliver the declared density part.|Potential PLY are meaningful but the common name is not a guarantee of both potential and density.|
| **XES** |RXES goes into the common absorption-rate branch; the output is marked XAS and uses the same name convention.|The result cannot be identified as a verified emission spectrum without reference physical comparison. It's not a setholder off.|
| **TDOBS** |Native path leads to generic DOS/observable loops; special time-related observational implementation has not been documented.|The name is not a guarantee of separate physical task. Cache fallback and subsequent rewrite flags must be understood by GEN.|
|**JXC: common FREL tensor with symmetrical equivalent neighbors; cluster SP-SREL path**|Rotation is disabled OR 1.EQ.1; cluster scratch files are not prepared in normal flow.|A simple periodic case can work, a general tensor/cluster workflow needs repair.|
|**SOCPAR: exchange output**|Unable to convert Ry→eV and test sign ABS(... )<0.|Values with eV label cannot be directly used; functional SOC part does not automatically disappear.|
|**PSHIFT: resonance fracturing**|The sum is used as well as the unearthed second channel and the interpolation does not control the bordered transition pi/2.|Phase shifts make sense, the derived splitting values are generally not reliable.|
|**PLOTPRS /TUTORIAL CHR-POT: spin and orbital curve resolution**|Both curves are entered into the same dataset; limits do not take into account the orbital curve.|Potential and charge graphs are applicable; faulty magnetic graph is the problem of Writer.|
|**MECHECK: promised full set of checks and summary**|Another test is compile-time disabled, NMEFORM is overwritten and the cumulative loop is not performed.|Sub-ADA/NAB control works; you will not switch on the full declared diagnostic by input.|
|**MEPLOT: current photon and electron scan**|The curves are rewritten through the wrong dimension.|This combination gives incomplete output; the new definition rejects it. Doubt about SEFIN is a particularly unverified interpretation.|
| **APS: broadening** |Expanded networks are not adequately controlled against capacities and the width depends on absolute energy zero.|The standard Lorentz breaking is not verified; without Broadening it can be basic APS applicable.|
|**POSANIPREP /POSANI: general data path NQ≠NT**|Writer/reader uses the number of types instead of the number of places.|Missing or outside the range of read data will not solve the general input validator.|
|**CHILANDAU: MPI and RDTAUFIL use**|Allocation/reduction errors and unprepared energy content.|Serial calculation without this reuse path is a different, more useful branch.|
|**COMPTON: EMD2D with non-zero PNMAX; some even shifted meshes**|Double scaling shifts the cut to another location; the construction loses the extreme points.|The result does not correspond to the declared cut/meshe. The definition protects entry restrictions.|
|**BAND: default EIGA/B**|Safe interruption is missing and continued until possible LOG(0).|Successful eigenvalue branch can work; solver failure must treat native code.|
|**CLPES: NE=1 with multiple types**|ECOR may remain from the last type instead of the selected type.|In this case, the energy axis may be incorrect; general solid geometry is not a problem itself.|
| **TUTORIAL: matching RWF-L** |Matching curve uses blank data.|Radial solutions and matching diagnostics are not the same; the graph needs repair.|
|**TUTORIAL: CORE-ALL outside the first suitable type/mesh/core case**|Incorrect JRWS index and unprepared variables in missing core states.|It is not generally completed © all core states of all types.|
|**TUTORIAL: POT-DECOMP/POT-XC for mixed magnetic types**|NONMAG is not renewed between types.|The resolution of output may be incorrect; comparison XC on fixed density is otherwise a deliberate example.|

Unprotected zero normalization FMAG/GILBERT/POSANITAU and cache or capacity errors AES are not proof that the entire main solver is not ready. These are specific faulty cases: there is a meaningful way in which conditions are met, the result cannot be trusted in their breach. The range is given in 2A and in the following error register.

Section 3 leaves separately tasks with directly documented breakage of the normal main workflow and overview of defective cases. The categories 2B and 3 are therefore not two competing conclusions: 2B specifies the applicability of the former caveats, sections 3 presents specific evidence and necessary corrections.

## 3. These tasks or specific branches are broken

### The intended normal workflow as a whole

| Task |Imported Problem|Evidence in the Source|Necessary correction|
|---|---|---|---|
| **FSCAT** |The error index of the harmonics with excess +L may exceed the field; the phase is calculated before treatment with a zero amplitude.| fscat.f:131,207–209; calc_chplm.f:155–162 |Fix mapping LM, filter inadmissible/zero CG folders and protect the phase from 0/0.|
| **T1** |Signed kappas are used as positive numbers m_j; negative channels fall out and indexing the whole base is broken down.| relt1.f:215–245 |Fix kappa degeneration/index and verify relaxation speeds. PRINT/RUNELOOP workaround only solves other initialization errors.|
| **NRAES** |IAES is zeroed by a shorter range than the allocated field; the convolution is added to the unitialized areas of different channels.| nraes.f:99,814 |Zero whole (2*NEMAX-1)*[NLMAX+1)*3 and add checks DOS metadata.|
| **BSDOSIJ** |J_ITJT is stored with energy indices into pair dimensions and read in a different order. It may cross the field and issue incorrect data.| bsdosij.f:319–349,434–436 |Unite dimensions/indexs, component offsets, rotations and final spectral output.|
| **BONDING** |CLUBONDING reads scratch 18/19 which does not prepare normal TAU_DRIVE with KTAUIJ=0.| tau_drive.f:70; tau_real_space.f:657–674,741–766; clubonding.f:78–108 |Connect actual exports TAUIJ to the consumer; then fix the centre, dimensions and output.|
| **NUCSSCPL** | Arrays are deallocated before their first allocation; the whole ME array is zeroed again for every type, and the mesh index is overridden to 1. | nucsscpl.f:115–133,341,366–368 | Fix array lifetimes, type-local initialization, and the mesh; NE=1 is not a fix for an integration task. |
| **XRS** |SIG(NEMAX+50,...) is zeroed only as NEMAX; the later sum uses the full range.| xrayspec.f:636–639 |Zero the entire field, correct the additional mesh and check the third integration capacity.|
| **BIS** |INITVAR recognizes it and GEN calls VBPES, but the routine requires TASK(3:5)=PES that BIS does not meet.| initvar.f:264; gen.f:1445; vbpes.f:174 |Implement the corresponding inverse-photoemission path. Alias BIS→VBPES would only mask the error and change meaning.|
| **SIGMA-XAS** |The calculation branch exists but the standard parser does not set internal TASK=SIGMA-XAS. Keyword SIGMA is not prefix recognition.|rsp.f:941; sigme.f:499; initvar.f and SECTION_FIND_KEYWORD|Complete the actual input path and verify XMCD workflow. Adding SIGMA switches to conductivity, not this branch.|
|**PLOTRHO separately**|RHORSQ is allocated only for the potential branch, but the density writer uses it without it.| fpplotprs.f:62,225,254,280 |Allocate RHORSQ independently in density branch. PLOTPRS with potential, this particular problem goes around.|
| **PLOTJALF / PLOTJORB** |FULLPOT ends up in an experimental PLY script before reading the current file; ASAwriter current ignores.| fpplot2dpr.f:71,194; plot2dpr_fp_ply.f:926; plot2dpr_asa.f |Fix dispatch/end and actual current reader. These are CONTROL selectors, not separate new implementations TASK.|

### Broken cases within otherwise meaningful tasks

| Rodina |Specific defective branch|
|---|---|
| PLOT2DRHO |ASA provides potential, not density. In FULLPOT the declared selectivity does not work and the common grid writer is not achieved.|
|POSANIPREP / POSANI data path|NQ≠NTwriter/reader may skip site or cross the field. The valid input will not correct it.|
| CHILANDAU |MPI and audited RDTAUFIL use path.|
| GILBERT, FMAG, POSANITAU |Relevant zero denominators; change ImE or number of points is not a general correction of normalization.|
| PSHIFT |Resonance splitting for missing second channel or without limited pi/2.|
| SOCPAR |Exchange column interpreted as eV without necessary transfer.|
| COMPTON |EMD2D with non-zero PNMAX and some even shifted meshes; rejects the definition of defective combination.|
| TUTORIAL |Matching plot RWF-L; CORE-ALL without core states / with non-compliant mesh; sticky NONMAG in POT-XC/POT-DECOMP for mixed types.|
| BAND |Failure of eigenvalue solver without safe continuation; possible LOG(0).|

**Correcting the conclusion of the previous selection audit:** audit PLOT-selectors provisionally designated PLOTRHO as fully meaningful. This has overlooked the previously detected missing allocation RHORSQ. The recontrol of the source confirms it; the separate PLOTRHO is correctly broken here. Its four definitions tests and no successful serialization are not a test of this native writer.

## 4. What changed in ASE2SPRKKR and what remains unresolved

### Input definitions and native-code fixes are separate layers

The earlier steps have added new definition, corrections to existing definitions, validators, descriptions and regression tests. These patches haven't been repairing Fortran***. U broken tasks thus the definition InputParameters can be syntactically correct, while the solve remains unusable for the declared purpose.

The existing definitions have been corrected for example for SCF, DOS, BSF, GILBERT, SIGMA, TORQUE, XAS, PHAGEN, JXC and ARPES; CHI and CHIXAS have received the audited new definition. Separate modules or shared families were added for further audited tasks. TUTORIAL has one complete definition of all 14 branches; older partial steps replace the final complete patch. TDOBS has a definition, but its native path remains generic DOS/observable loop: no specific time-based observational calculation was documented. They are therefore among the **working with a large ALE in the name and expected meaning**.

Audit-only conclusions BIS, SIGMA-XAS, VBUPS/VBXPS and PLOTJALF/PLOTJORB do not provide fake new solver definitions: the first two branches are not functionally achievable by common input, the aliases already share VBPES and current graphs block native dispatch.

### The ARPES family

ARPES, AIPES, SPLEED, BAND and CLPES use KKRSPEC but have different IBLOCH paths, models and parametric meanings. VBPES/VBUPS/VBXPS, by contrast, useGEN/VBPES. They are not two interchangeable implementations of the same complete calculation.

Delivered warning for old TASK=AIPES/BAND/SPLEED within definition ARPES; started at validation/saving. It's DepressionWarning, standard Python filters can hide it. Original ARPES does not deprecate itself. Common KKRSPEC basis instead of import of whole ARPES definition is suitable next refactor, not yet completed.

The difference remains: new specialised definitions use correct **CRYS_VECS***, original ARPES still offers incorrect **CRYS_VEC***. Plural native keyword also switches the use of crystal primitive vectors according to its presence; it is not a common positive boolean within the meaning of the name. Native IMV_INI v Ry erroneously writes VIH, while eV variant sets correctly VIL. These facts cannot be corrected by depreciation warning.

### TUTORIAL as a complete task

Final definition includes: **YLM, JLX, JLX-ASYMPT, CHR-POT, POT-DECOMP, POT-XC, RWF-L, RWF-E, CORE-ALL, DOS-FEG, DOS-SS, DOS-CLU, DOS-BULK, GEN**. This is not one physical calculation with 14 names. Defaults and supported scheme are selected by branch: in particular RWF uses scalar relativity and actual radial mesh, DOS-CLU must have cluster and POT-XC compares fixed density functionals.

## 5. Defaults: explicit input values, omission, and the actual fallback

**Rule of this audit:** default which is omitted from the input and represents a known value must correspond to an effective native fallback. `None` may indicate a delegation to EF, potential or metadata; it is not a claim that native value is Python None. The value listed may be different, but the difference must have the right meaning. Parameter which the solver later overwrites is not an effective user control only because it is registered.

### Important shared findings

|Case|Correct Meaning/Reparation|
|---|---|
| EMAX sentinel |INIT_MOD_ENERGY converts the value less than -99 to EF. **EMAX=-1 is not EF* but absolute -1 Ry. Deleted EMAX after the relevant fallback may mean actual EF.|
|Relative vs. absolute energy|General EMINEV/EMAXEV are eV offsets against EF. EKREL TASK bounds are absolute Ry. SPLEED incident EV is absolute incident energy; does not transfer here ARPES meaning.|
| Presence-sensitive options |NSHLCLU=0 present is not the same as omitted NSHLCLU present. NEPHOT=1 present is not the same as single-photon fallback without keyword.|
|GRID and NE|GRID=2 s NE=1 may divide by zero. GRID=3 s NE=1 has a usable first point but uninitialized integration weights; it is suitable only where it is not integrated. GRID=0 has its own single-energy meaning.|
| ENERGY vs. TASK |COMPTON needs SPLITSS in both sections. EKREL search NE/EMIN/EMAX are in TASK. OPTABS actual own network controls TASK parameters.|
| Potential context |NL, IT range, existence shells, NQ/NT, KMROT, moments and actual FULLPOT format not generally detectable from InputParameters without POTFIL.|

### Explicitly written differences that are not errors in themselves

| Task |ASE2SPRKKR / patch vs. native behaviour|Assessment native default|
|---|---|---|
| SCF |NE=32 versus native 30; some custom mixing options. Omitted MIXOP remains natively 0.2 even if another MIX is read later.|30 can be cheaper, not automatically converted. MIXOP cannot be described quietly as current MIX.|
| DOS |Explicit GRID=3, NE=300, EMIN=-0.2, EMAX=1 Ry. Deleted native upper limit would be EF.|EMAX=1 is a meaningful scan of unfilled states; EF would limit output to occupied.|
| BSF EK |Explicit EMAX=-1 Ry has been preserved.|Usually better is EF or your own growing interval. It shall not be marked -1 as a sentinel EF; at EMIN=-0.2 it shall be a descending scan.|
| BSF |NK=300 against native 800; KK EMIN=0.7 Ry.|The density should be converged; KK default is not guaranteed by Fermi cut.|
| XAS |Explicit EMIN=-0.2 Ry against native EF; EMAX=4 has been corrected to actual member.|EF is more natural for standard absorption, but native XAS may still overwrite EMIN. Check EMF and effective limits.|
| PHAGEN |NE=32 versus 30, ImE=0 versus 0.01; NITER/MIX are converted to 1/1 when exported.|Do not include in the description the number of iterations from registered common default. Convergence of contours has not been compared.|
| JXC / FMAG |The preserved explicit ImE=0, even if the generic native EIMAG is 0.01.|U GRID=5 decides the whole complex conture; the zero itself does not prove usability or better choice.|
| DMI / CHILANDAU / CHI |GRID=5, NE=30, ImE=0.01 taken over by native integration workflow. CLURAD=2.2 is an explicit practical option.|Native number of points is appropriate default base; radius is not a guarantee of convergence.|
| SIGMA | ImE=0.01 proti surface native nule. |Smearing significantly affects transport; none of the values is universally better without convergence of study.|
| GILBERT |Fixed GRID=0/NE=1 for EF, non-efficient bounds removed; NL=4 is explicitly written.|Native One-energy path is better than an inadequate linear contour.|
| OPTABS |TEMPEL=300 K/NKT=6 versus reader values 1e6/99; NEB=30/NOMEGA=100, OMEGA1=0/OMEGA2=8 eV.|300 K is a more practical physical start. Native extreme temperature is not a reason for automatic take-over; validate actual temperature use.|
| ARPES |NE=300, EV=-8..5, EWORK=4.2 eV, damping 0.05/2 eV, NL=4 and layers 50/57; native has other numbers/damping.|This is a written option, not a omitted bug. Unable to declare native 200 points and other layers as better without convergence.|
| AIPES / SPLEED / BAND |AIPES valence EMAXEV=0; SPLEED owns 11..60 eV/200 points and 0.2/0.2 eV; BAND explicit zerodamping and EIG_IME.|Correct task-specific differences; do not unite blindly according to ARPES.|
| KKRSPEC N_LAYDBL |Explicit two components bypassing native typo double initialization LANZ1 instead of LANZ2.|Absent ***default is not safely initialized; explicit values are better for repair native.|
| T1 |PRINT=1 and RUNELOOP=true are explicit workarounds of defective initialization.|Native defaults are not better here. Ani workaround does not fix kappa sum. |
|Diagnostic fences|GRID/NE often only provide a common startup allocation; the last selectors use GRID=0/NE=1/ImE=0.|More points do not improve the radial or spatial graph. Do not attribute startup energy to physical spectral meaning.|
| APS/NRAES/AES/VBPES |Explicit NE=300 is the capacity or assumed mesh; actual NE is taken from DOS/TAU.|Better future solution is deriving and validation from the real file header including reserve breaking.|
| WFPLOT / MEPLOT / MECHECK |Single-point diagnostics written; context EMAX and NL are omitted. WFPLOT NL=3 with writes, heavy context may have 4.|Higher NE/NL change by purpose; do not make one example generic native default.|

Final TUTORIAL definition writes different energy defaults: mathematical branches GRID0/NE1/ImE0; RWF GRID3/NE10/ImE0 with absolute 1e-7..1 Ry; DOS/GEN GRID3/NE300/ImE0.01 s -0.2..1 Ry. It's an explicit practical election. STANDARD and relevant BZINT are also recorded by branch.

## 6. MPI and outputs

`mpi=True` in the definition means declaration support, not performed MPI regression run `mpi=False` may be a consequence of documented native guard or conservative limitation of unverified workflow; these reasons should not be mixed.

| Skupina |Result of source audit|
|---|---|
|SCF, standard DOS, BSF, EKREL, JXC, DMI, COMPTON, POSANI, FMAG|There is a source justification MPI for the audited path; numerical serial/MPI match not measured. FMAG completes DENSITY before general GEN guard.|
|XAS/XES/XRS, STONER, VBPES, APS/NRAES/AES, JALF/JORB, diagnostic GEN fences|A general energy loop is not a guarantee of completed postprocessing. Relevant native guards or missing reductions justify serial interfaces.|
| CHILANDAU, NUCSSCPL, PHAGEN |Specific source errors in MPI/locations/barriers.|
| GILBERT |The existing definition has an implicit MPI=true; it cannot be claimed that the patch supplied generally set mpi=False.|
|SIGMA, POSANIPREP, POSANITAU, new KKRSPEC wrappers|Serial option is in part conservative audit boundary, not proof that there is no MPI in the source. BAND owns MPI implementation actually has.|

U many new definition are not specialised Python output parser. Common output reader is used. Grace/PLY or rate file is therefore not automatically structured physical result ASE2SPRKKR.

Fixed names `ME`, `MEi`, `sigmaXXXXX.dat`, `curdns1..3`, `fort.*`, some PLY and `arpes.dat` are not always isolated using CONTROL.DATASET. U AES/PHAGEN/VBPES/POSANIPREP additionally existing files can change the cache path. Numerical heads are not always free ENERGY.NE: AES/NRAES have 2*NE-1 spectral points, XRS adds additional mesh and JORB/JALF write radial/angle data.

## 7. Tests, revisions, and verification limits

The tests verified the discovery definitions, text and file roundtrip, copy/pickle, presence or omission defaults and rejection of specific invalid combinations. They did not allow numerical SPR-KKR and did not verify spectra, convergence, cache origin or physical reference results.

Most of the separate patches were formed against ASE2SPRKKR `refactor` **267512c4b8ccbbfd37af7ede898febc066b7e3ee**. The oldest WFPLOT base used **8549af97ae106b08749ac1f2dc4b819a72b988e1**. It is not claimed that these changes are located on today's head of external repository.

**Patche were mostly independent. Counts of the whole set are not counted.** TUTORIAL complete replaces partial 0049–0052. ARPES audit and warning have separate patches; their common versions have also been tested. All task patches together were not integrated or tested together. The consolidation of the audit does not implement this integration.

Typical command running in the appropriate checkout with the Python test:

```bash
PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 python -m pytest src/ase2sprkkr/input_parameters/test -q
```

For today's purely documentary consolidation, the solver has not been restarted or made new changes to the definitions. A new pytest that would only test the content of this table would not add verification native errors.

### Recorded results of the individual deliveries

Below are historical test results corresponding to the definition above its base, not a new summary set. The test names are the names of modules delivered; the results are from individual audits. All canonical individual audits are part of the support package in order to remain traceable for their specific scope.

|Delivery|Test module / range|Whole input set: passed|Detailed background in ZIP|
|---|---|---|---|
| WFPLOT | `test_wfplot.py` | 128 | [WFPLOT-audit.md](audit-evidence/WFPLOT-audit.md) |
| MEPLOT | `test_meplot.py` | 128 | [MEPLOT-audit.md](audit-evidence/MEPLOT-audit.md) |
| MECHECK | `test_mecheck.py` | 126 | [MECHECK-audit.md](audit-evidence/MECHECK-audit.md) |
| SOCPAR | `test_socpar.py` | 125 | [SOCPAR-audit.md](audit-evidence/SOCPAR-audit.md) |
| PSHIFT | `test_pshift.py` | 123 | [PSHIFT-audit.md](audit-evidence/PSHIFT-audit.md) |
| FSCAT | `test_fscat.py` | 121 | [FSCAT-audit.md](audit-evidence/FSCAT-audit.md) |
| PLOTPRS | `test_plotprs.py` | 120 | [PLOTPRS-audit.md](audit-evidence/PLOTPRS-audit.md) |
| PLOT2DPR | `test_plot2dpr.py` | 124 | [PLOT2DPR-audit.md](audit-evidence/PLOT2DPR-audit.md) |
| T1 | `test_t1.py` | 118 | [T1-audit.md](audit-evidence/T1-audit.md) |
| EKREL | `test_ekrel.py` | 117 | [EKREL-audit.md](audit-evidence/EKREL-audit.md) |
| CLXPS | `test_clxps.py` | 121 | [CLXPS-audit.md](audit-evidence/CLXPS-audit.md) |
| STONER | `test_stoner.py` | 120 | [STONER-audit.md](audit-evidence/STONER-audit.md) |
| VBPES | `test_vbpes.py` | 120 | [VBPES-audit.md](audit-evidence/VBPES-audit.md) |
| APS | `test_aps.py` |116 (original patch); 9 isolated tests after repair n≤5| [APS-audit.md](audit-evidence/APS-audit.md) |
| NRAES | `test_nraes.py` | 116 | [NRAES-audit.md](audit-evidence/NRAES-audit.md) |
| AES | `test_aes.py` | 121 | [AES-audit.md](audit-evidence/AES-audit.md) |
| FMAG | `test_fmag.py` | 116 | [FMAG-audit.md](audit-evidence/FMAG-audit.md) |
| JORB | `test_jorb.py` | 117 | [JORB-audit.md](audit-evidence/JORB-audit.md) |
| JALF | `test_jalf.py` | 118 | [JALF-audit.md](audit-evidence/JALF-audit.md) |
| BSDOSIJ | `test_bsdosij.py` | 121 | [BSDOSIJ-audit.md](audit-evidence/BSDOSIJ-audit.md) |
| BONDING | `test_bonding.py` | 120 | [BONDING-audit.md](audit-evidence/BONDING-audit.md) |
| NUCSSCPL | `test_nucsscpl.py` | 122 | [NUCSSCPL-audit.md](audit-evidence/NUCSSCPL-audit.md) |
| TORQUE | `test_torque.py` | 119 | [TORQUE-audit.md](audit-evidence/TORQUE-audit.md) |
| XES | `test_xes.py` | 124 | [XES-audit.md](audit-evidence/XES-audit.md) |
| XAS | `test_xas.py` | 123 | [XAS-audit.md](audit-evidence/XAS-audit.md) |
| XRS | `test_xrs.py` | 124 | [XRS-audit.md](audit-evidence/XRS-audit.md) |
| PHAGEN | `test_phagen.py` | 121 | [PHAGEN-audit.md](audit-evidence/PHAGEN-audit.md) |
| JXC | `test_jxc.py` | 122 | [JXC-audit.md](audit-evidence/JXC-audit.md) |
| DMI | `test_dmi.py` | 122 | [DMI-audit.md](audit-evidence/DMI-audit.md) |
| COMPTON | `test_compton.py` | 130 | [COMPTON-audit.md](audit-evidence/COMPTON-audit.md) |
| POSANI | `test_posani.py` | 132 | [POSANI-audit.md](audit-evidence/POSANI-audit.md) |
| POSANIPREP | `test_posaniprep.py` | 121 | [POSANIPREP-audit.md](audit-evidence/POSANIPREP-audit.md) |
| POSANITAU | `test_posanitau.py` | 123 | [POSANITAU-audit.md](audit-evidence/POSANITAU-audit.md) |
| DOS | `test_dos.py` | 120 | [DOS-audit.md](audit-evidence/DOS-audit.md) |
| GILBERT | `test_gilbert.py` | 122 | [GILBERT-audit.md](audit-evidence/GILBERT-audit.md) |
| CHILANDAU | `test_chilandau.py` | 120 | [CHILANDAU-audit.md](audit-evidence/CHILANDAU-audit.md) |
| CHI | `test_chi.py` | 119 | [CHI-audit.md](audit-evidence/CHI-audit.md) |
| OPTABS | `test_optabs.py` | 121 | [OPTABS-audit.md](audit-evidence/OPTABS-audit.md) |
| SIGMA | `test_sigma.py` | 122 | [SIGMA-audit.md](audit-evidence/SIGMA-audit.md) |
| CHIXAS | `test_chixas.py` | 121 | [CHIXAS-audit.md](audit-evidence/CHIXAS-audit.md) |
| BSF | `test_bsf.py` | 119 | [BSF-audit.md](audit-evidence/BSF-audit.md) |
| SCF | `test_scf.py` | 123 | [SCF-audit.md](audit-evidence/SCF-audit.md) |
| ARPES | `test_arpes.py` | 122 | [ARPES-audit.md](audit-evidence/ARPES-audit.md) |
| SPLEED | `test_spleed.py` | 119 | [SPLEED-audit.md](audit-evidence/SPLEED-audit.md) |
| AIPES | `test_aipes.py` | 119 | [AIPES-audit.md](audit-evidence/AIPES-audit.md) |
| BAND | `test_band.py` | 118 | [BAND-audit.md](audit-evidence/BAND-audit.md) |
| CLPES | `test_clpes.py` | 116 | [CLPES-audit.md](audit-evidence/CLPES-audit.md) |
| TUTORIAL complete |`test_tutorial.py (final complete version)`| 148 | [TUTORIAL-complete-audit.md](audit-evidence/TUTORIAL-complete-audit.md) |
| BIS |`No new definition; source audit`| — | [BIS-audit.md](audit-evidence/BIS-audit.md) |
| TDOBS | `test_tdobs.py` | 113 | [TDOBS-audit.md](audit-evidence/TDOBS-audit.md) |
| SIGMA-XAS |`No new definition; source audit`| — | [SIGMA-XAS-audit.md](audit-evidence/SIGMA-XAS-audit.md) |
| PLOT-selectors | `test_plot_selectors.py` | 112 | [PLOT-selectors-audit.md](audit-evidence/PLOT-selectors-audit.md) |
| PLOT2D-selectors | `test_plot2d_selectors.py` | 115 | [PLOT2D-selectors-audit.md](audit-evidence/PLOT2D-selectors-audit.md) |
| VBUPS-VBXPS |`No new definition; source audit`| — | [VBUPS-VBXPS-audit.md](audit-evidence/VBUPS-VBXPS-audit.md) |
| PLOTJALF-PLOTJORB |`No new definition; source audit`| — | [PLOTJALF-PLOTJORB-audit.md](audit-evidence/PLOTJALF-PLOTJORB-audit.md) |
| ARPES warning |Single warning set and together with ARPES patch|112 separately; 126 together|Delivery 0047; standard ARPES is not monitored|

## 8. Priorities for further fixes

1. **SPR-KKR before physical use:** correct fields/indexs FSCAT, NRAES, BSDOSIJ, NUCSSCPL and XRS; kappa sum T1; bitrate BONDING; writer PLOTRHO. These are documented flaws that the Python definition will not solve.
2. **Reachability and significance:** resolve BIS and SIGMA-XAS; separate selective PLOT2D/current writer from experimental PLY termination; physically compare XES and specify AIPES model.
3. **ASE2SPRKKR:** correct original ARPES CRYS_VEC, unify actual KKRSPEC base and distinguish serial metadata from native ban. Judge BSF default EMAX=-1 and confirm the change in the energy window as a change in public behaviour.
4. **Possibility/cache:** validate TAU/DOScapacity, actual NQ/NT, meshes, shells and magnetic geometry where the file is available. Fixed NE or positive IT do not replace this check.
5. **Integration:** join a separate patch on one familiar revision, solve shared helpers and overlays and run one common input_parameters set. Only then deliver representative numerical and serial/MPI regression cases for supported physical pathways.

## 9. Disabled tasks — the combined list at the end

Only the list of branches from the submitted revision; no other pretend definitions were created for them. U public blockage is explicitly supplied PUBLIC_VERSION. U NEGF is the off old path over KKRRSP, not by this conclusion the whole separate program KKRNEGF.

|Task / Branch|Reason|Place|
|---|---|---|
| **FORCETEN** |UNDER_CONSTRUCTION before calculation| gen.f:1055–1057 |
| **XMO** |Public version: not available / UNDER_CONSTRUCTION| xrayspec.f:283–287 |
| **OPTICS** |Public version: not yet available| rsp.f:887–888 |
| **MAGNET** |LINRESP_MAGNET is a tube with STOP| linresp_magnet.f:1–4 |
| **CHIDYN** |LINRESP_SUSC_DYN is a tube with STOP| linresp_susc_dyn.f:1–4 |
| **PHONONS** |Public version explicitly stops task| check_input.f:34–36; rsp.f:962–963 |
| **PUMP** |LINRESP_PUMP is a tube with STOP| linresp_pump.f:1–4 |
| **OPM** | OPM scheme not available | check_input.f:136 |
| **SIGMA-DNS** |SIG_DNS to UNDER_CONSTRUCTION| sig_dns.f:13 |
|**NEGF over KKRRSP**|From version 8.3 refers to separate KKRNEGF and stops this path| rsp.f:905–907 |
