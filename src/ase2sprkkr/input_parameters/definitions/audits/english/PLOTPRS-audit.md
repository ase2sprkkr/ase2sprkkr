# PLOTPRS — SPR-KKR 9.7 and ASE2SPRKKR audit

PLOTPRS and its selective plotting variants are **prepared in `tasks-workflow`**, not available in `develop`. The native SPR-KKR 9.7 plot routines have separately documented limitations.

## Interface and range

`initvar.f:187–190` at keyword PLOTPRS turns on PLOTPRS(0:3), i.e. all three graph groups. `gen.f:1460–1463` calls FPPLOTPRS. The actual routine `fpplotprs.f` does not read task-specific parameters. The definition therefore does not add choice of atom, orbital, radial mesh, interval r or group switches that the source would not read.

Use: `InputParameters.create('PLOTPRS')`. The definition uses CONTROL, TAU, MODE, STRCONST, ENERGY, TASK, SITES; executable=kkrgen, mpi=False. Routine writes common files without distribution of work between wounds; MPI nor numerical startup has been verified here.

| Grafy | ASA | FULLPOT |
|---|---|---|
|Potential|spherical r*VT and BT for each atomic type|in addition VNST and BNST by non-spherical components|
| Hustoty |Spherical charge density, spin/orbital density under magnetic mode|additionally non-spherical RHO2NS|
| Shape functions |No output SFN|FLMSF for each radial mesh IM=1..NM|

The spherical charge curve is RHOCHR*r2/(4*pi), not free RHOCHR (`fpplotprs.f:254`). Radius is in a.u.; potential labels distinguish r*V and B. Shape functions are coefficients of existing radial representation, not 3D visualization of cell shape. NONMAG is taken from a common CONTROL and limits the number of panels.

The purpose of task is to draw data read into radial fields, not SCF to calculate new density. FPPLOTPRS samo shape functions does not create. Initialization FULLPOT has a separate way of reading/creating SFN; the availability of correct representation depends on the potential, mesh and file SFN, not on the new parameter of this routine. PLOTPRS works also in ASA where the condition `PLOTPRS(3).AND.FULLPOT` shape-function graphs skip (`:392`). There is therefore no reason to introduce a general ban ASA.

## Defaults and Initialization

| Option |Definitions|Justification|
|---|---|---|
| TASK |explicit PLOTPRS|natively turns on all three groups|
| ENERGY.GRID |explicit 1|simple initialization path; does not specify radial graphs|
| ENERGY.NE |explicit 1|INIT_MOD_ENERGY requires NEMAX>0 even when energy loop is switched off|
| ENERGY.ImE | 0 |not used as smearing radial curves|
| ENERGY.EMIN |explicit -0.2 Ry|identical s native EMINDEFAULT; for inefficient graphs|
| ENERGY.EMAX | None, is omitted from the input |native fallback to EF|
| SITES.NL | None, is omitted from the input |context cutoff from potential|
| CONTROL.SFNFIL |Common default DATASET.sfn|can explicitly enter your own filename|

`initvar.f:211–218` turns off INITELOOP and RUNELOOP. Yet the main KKRGEN calls INIT_MOD_ENERGY; it sets NEMAX from NETAB and ends the programme at zero (`init_mod_energy.f:192–197`). It would therefore not be safe to completely omit ENERGY. NE=1 is a written technical option ASE2SPRKKR, not a statement of native omitted default. The higher number of points does not add radial graph points; the definition draws attention to it. Default NE=300 used in other graphs would not be better here.

Unwritten EMAX and NL maintain native context. CONTROL, TAU and STRCONST are shared definitions; task does not evaluate scattering-path operator when graphing. Explicit general options cannot therefore be interpreted as additional parameters FPPLOTPRS.

## Restrictions and selective variants

GEN forbids DMFT and LDAU and for PLOTPRS (`gen.f:211–222`); the definition rejects them. In the source there is no general prohibition FULLPOT nor requirement IREL=3, therefore the definition allows FREL, NREL, SREL and SP-SREL. The numerical combination of all modes has not been verified.

TASK.HFF may later in INITVAR set IHYPER but for PLOTPRS does not run the usual energy loop nor calculate hyperfine results. The definition draws attention to an ineffective requirement; it does not claim that the flag is not loaded at all.

SPR-KKR recognises also PLOTPOT, PLOTRHO and PLOTSFN (`initvar.f:192–196`). These selective variants ** are not registered with this patch**: task is PLOTPRS with all groups. It is an explicit limit of support, not an undisclosed equivalence of all variants. In addition, the separate PLOTRHO is incorrect in the source, see below. PLOTPRS cannot be forced to shut down groups by combining fake flags: its keyword turns them all on.

## Errors and risks found SPR-KKR

1. **Standalone PLOTRHO uses a non-loqued field.** RHORSQ is allocated only inside PLOTPRS(1), when drawing potential (`fpplotprs.f:62`). The density branch allocates only Y (`:225`), but writes RHORSQ (`:254`) and transfers it XMGRTABLE (`:280`). If only density is switched on, RHORSQ does not have an allocation. Normal PLOTPRS passes this error because the first branch also runs; does not mean means that the selective option is correct.
2. **Spin and orbital density use the same chart dataset.** Both calls XMGRTABLE have panel=1 and curve index=0 (`:282–287`). So the format gets two tables with the same goal instead of separate curves. The specific import effect depends on the display; the data are not clearly marked as two different curves.
3. **Limities of the second density panel ignore orbital density.** Y1MIN/Y1MAX are determined only from RHOSPN (`:257–258`) although with writes also RHOORB. The orbital curve may extend beyond the set range.
4. **Great L does not fit into format in some non-spherical legends.** Potential and density branches use `99001 FORMAT(I1,',',I2)` for l,m. For l>=10 there is an asterisk instead of l; SFN the branch, on the other hand, selects a wider format according to size l,m. A general ban on high cutoff due to a cosmetic error has not been introduced.

RHORSQ is local allocatable; missing explicit DEALLOCATE at the end is not in itself proof of memory leakage because Fortran releases local allocatables on return. The audit therefore does not call this detail another mistake.

The output uses DATSET0, derived from CONTROL.DATASET in the main programme (`kkrgen.f:239`). Files are created over XMGRHEAD as V, VNS, RHO, RNS and SFN_M by group. Patch does not contain specialized output parser graphs.

## Verification scope

The prepared Python definition has test coverage for serialization and selected validation rules. Counts from old standalone patch runs are not the test results of the current `tasks-workflow` branch. No numerical native SPR-KKR reference calculation is claimed here.
