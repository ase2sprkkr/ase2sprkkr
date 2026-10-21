# COMPTON — current definition audit

`compton.py` implements the input contract of periodic FREL Compton calculations using two GRID=5 contours and SPLITSS. The native branch is active but has documented limitations; no numerical validation is claimed.

## State and scope

The current definition supports periodic FREL Compton momentum densities through `kkrgen` on the documented split-contour path. Native alternatives are assessed separately below.

## Outputs and filenames verified from SPR-KKR

Native `compton.f:1314–1354` writes `DATSET_BND.mcp`, `DATSET_COR.mcp`, and `DATSET_TOT.mcp`, plus `DATSET_data.mcp` for `PRINT>0`. `MD2D` selects `_mcp2d.dat`, whereas `EMD2D` selects `_emd2d.dat`; both write the `BND`, `COR`, and `TOT` components (`compton.f:1439–1458`). These represent band, core, and total contributions. `MD2D` is an integrated projection; `EMD2D` a 2D momentum-density cut.

`kkrmain.f:239–245` first sets `DATSET0=DATSET`, then appends nonempty `ADSI` to `DATSET`. COMPTON results consequently include ADSI in their names. Current `compton.py::output_files()` uses `CONTROL.DATASET()` alone, and thus declares the wrong paths for nonempty ADSI. This is a manifest discrepancy, not a native calculation defect.

Native `init_mod_compton_pmesh.f:183–185,269–279` scales `PNVEC` first by `DPN`, then again by `PNMAX` for `EMD2D`. This confirms the documented limitation of nonzero longitudinal cuts.

## Source connections and inputs

- `initvar.f:480–489`: COMPTON switches RDTAUDIR to true and normal energy
loop to false. `gen.f:1572–1577` calls COMPTON.
- `compton.f:167–170`: requires IREL=3 and compatible constant NSYMMAX.
- `kkrgen.f:264–274`: default TAUDIR is DATASET.tau.dir; can be entered TAUDIR
a TAUPATH in the TAU section. Both parameters are added to the definition.
- `init_mod_energy.f:284–302`: controls the number of energy and compatibility TAU
directories; some discrepancies switch on conversion, another number of energy ends STOP.
**Definition does not ensure creation or content of cache.** Calculation needs suitable
TAU backgrounds from the previous step. Python tests don't verify their compatibility.
- `compton.f:203–228`: USE_ZJOVZZ, SPLITSS, MEMODE and NPMAX.
- `init_mod_compton_pmesh.f`: PNVEC, PPVEC1, PPVEC2, PNMAX, PPMAX,
  NPN, NPP, MD2D, EMD2D, SHFT, SHFT_PP, SHFT_PN, INSIDE_SPHERE, USE_PEQUIV.
- `compton.f:263`: SQQ_RESOLVED. Supplemented as the choice of resolution by pair of places.

The definition deliberately covers the standard split-contour procedure. ONLY_CORE and USE_ZJOVZZ
they have other branches, they are not exposed or declared off. GRID=5/5 and POINTS
are limitations of the audited scope, not the claim that all other native
The energy grids are forbidden. The numbers NE shall be positive.

## Defaults: written to the input versus omitted

| Parameter |Definitions and writing|Native behaviour|Assessment|
|---|---|---|---|
| GRID, NE | writes 5/5, 30/50 |equal in INIT_MOD_ENERGY for COMPTON|Taken native defaults|
| ImE, EMIN | writes 0.01, -0.2 Ry |the same| native defaults |
| EMAX |not issued and is omitted from the input|set to EF for COMPTON|the same effective behaviour|
| ENERGY.SPLITSS |true, always writes keyword|Initialization for COMPTON sets true|Maintained|
| TASK.SPLITSS |true, always writes keyword|its absence in COMPTON overwrites SPLITSS to false|explicit option maintaining the two-konture procedure|
| MEMODE |1, writes| 1 |the same|
| NPMAX |250, writes|250; hodnoty <=0 returns to 250, >10000 cropped|Same default, validation 2..10000 prevents NP=1|
| PNMAX, PPMAX |0, 10 a.u., entered|the same|PNMAX=0 leads to a single longitudinal point, the validator warns|
| NPN, NPP |50, 50, entered|the same; for PNMAX=0 the NPN is effectively converted to 1|parameter and effective size vary|
| USE_PEQUIV |0, writes|0; outside 0..2 returns to 0|default same invalid values rejected|
|directional vectors|singles, omitted from the input|Normal z, transverse axis x/y|the same|
| MD2D, EMD2D, shift flags, INSIDE_SPHERE, SQQ_RESOLVED |false not to be published|false in the absence of keywords|the same|
| TAUDIR, TAUPATH |singles, omitted from the input|DATASET.tau.dir, current directory|the same|

PNMAX=0 is not reasonable default for the whole profile but can serve for a single point.
Native default is therefore left with notice; the user is to specify the range.
Native NPMAX=250 is a suitable start, needs convergence control; MEMODE=2
Allows to compare the calculation without interpolation. The default is not a guarantee of convergence.

## Mistakes and restrictions native code

### Two independent readings SPLITSS

INIT_MOD_ENERGY sets SPLITSS for COMPTON but `compton.f:210` again
Loads **from TASK*** where the absence of keyword means false. ENERGY.SPLITSS samo
not enough. Definition always writes both keywords and validator requires both true.
It is not automatic synchronization of two editable values: shutting down any
refuse to be validated. USE_ZJOVZZ additionally SPLITSS would force false (`:212`),
Therefore, this alternative is outside the scope of the definition.

### Directional vectors are non-ortogonalized

`init_mod_compton_pmesh.f:49–121` normalizes specified vectors and creates a third
cross product vector. It doesn't remove the longitudinal component of the second specified vector.
Zero or parallel vectors can divide by zero; non-circular creates a deformed base.
The definition controls the final nonzero vectors of the length 3 and the perpendicularness of the pair.
At PNVEC, PPVEC2 is ignored; the definition rejects this misleading combination.
Without PNVEC, both transverse vectors or none shall be entered.

### Even shifted grids lose extreme points

`init_mod_compton_pmesh.f:206–245` uses ±INT((N-1)/2) but original step
the grid is MAX/(N-1). For N=50 instead of 50 points is obtained 49 and the range is not reached
both intended regions. For N=51 the design is consistent.
The definition requires an odd number in the shifted direction and a maximum of one shift flag.
This is conservative protection of the range, not the claim that native even input always crashes.

### EMD2D: longitudinal scaling takes place twice

`init_mod_compton_pmesh.f:183–185` is already multiplying PNVEC by DPN.
In EMD2D then `:273–279` sets DPN=PNMAX and rescales the same vector.
The cut is not thus generally formed on the declared PNMAX. For example, for NPN=50,
PNMAX=2 without shift flags is the first range 2/49 and the resulting 4/49 instead of 2.
The definition allows EMD2D only for PNMAX=0 where this particular problem does not affect
the location of the cut. Other EMD2D entries are rejected. Native code has not been changed.

### Other significant ALE

- `NPMAX=1` results in a division NP-1 at interpolation (`compton.f:446`;
  `compton_mod1_me0.f:144`). The definition requires at least 2.
- `compton.f:185–194` at GF_CONV_RH switches the convention to ZJ and SOLVER_FP to BS.
I FULLPOT therefore changes the convention; immediate compatibility with cache should be checked.
- S INSIDE_SPHERE is PFEMAX=min(PNMAX,PPMAX); at native PNMAX=0 is cut
to 1e-10. This type of ball does not provide the usual profile.
- Specific systems require compatible sizes NKMQ; native controls
`compton.f:678` may end up STOP. The definition without potential will not guarantee this.

## MPI

`compton.f:581–592` requires NE for each active contour at least NPROCS.
Energy loop uses MPI_DISTRIBUTE and reduces results MPI_SUM_R
(`:1273–1296`). `mpi_distribute.f` divides energy backwards, so last
the point processes the rank 0; corresponds to the condition of the calculation of the core contribution (`compton.f:626`).
Therefore, mpi=True is the source justified. No actual MPI run was performed.
