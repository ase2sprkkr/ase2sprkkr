# AES audit and new definition

Base: currently fetched ASE2SPRKKR refactor 267512c4b8ccbbfd37af7ede898febc066b7e3ee and supplied sources SPR-KKR 9.7. V ASE2SPRKKR here AES does not have its own definition. The separate patch 0016-aes.patch contains aes.py and test_aes.py, does not include previous tasks.

## Workflow and input

Fully relativistic Auger electron spectroscopy in kkrgen. initvar.f:299–308 sets INITELOOP=true, RUNELOOP=false, RDTAUDIR=true, WRTAUDIR=false and USELPLS1=true. Routine AES reads initial TAU directory. This is not a postprocessing DOS file as NRAES.

TASK reads IT, CL, NEME, PRINTME and PRINTAME (aes.f:128–144). Options CALCLS, OUTGNU or GAMMA from nonrelativistic branches do not belong here. The source only works with NQ=1 (aes.f:230). The definition documents this restriction; the general validator InputParameters does not have potential metadata and therefore cannot verify the number of places. NQ=1 does not mean necessarily NT=1 for alloy. Explicit IT is to be an index of the type existing in potential.

check_input.f prohibits FULLPOT and requires IREL=3. GEN prohibits DMFT and LDA+U. The definition rejects incompatible MODE and both correlations; SPR-KKR must verify FULLPOT and NQ from the potential. MPI is disabled because AES is not included in the list tasks continued after MPI energy loop GEN.

## Defaulty

| Option |Native behaviour|Definitions| Serialization |
| --- | --- | --- | --- |
| TASK.IT | 1 | 1 | default omitted |
| TASK.CL | 2P | 2P | default omitted |
| TASK.NEME | 10 | 10 | default omitted |
| PRINTME / PRINTAME |false in absence| false | omitted |
| ENERGY.GRID / NE / ImE |without specific AES meshes; initial mesh must correspond TAU| 3 / 300 / 0,01 Ry |explicit|
| ENERGY.EMIN | -0,2 Ry | -0,2 Ry |explicit|
| ENERGY.EMAX | EF |None, delegation per source| omitted |

Omitted defaults TASK correspond to the source. GRID=3 and NE=300 are the explicit choice of the definition for initial mesh and capacity, not the claim of native AES defaults. When using an existing TAU, input mesh must correspond to the file; fixed 300 is not a detection of its content. A better solution would be to derive mesh and capacity from metadata specific TAU. Left out EMAX correctly leaves EF out of potential.

NEME is the gross net matrix elements, the actual number used is min(TASK.NEME, NE z TAU). The definition allows NEME>=3 because the interpolation uses adjacent points on both sides. ENERGY.NE>=3 only checks the allocation capacity, does not guarantee actual NE>=3 in the file.

## Source findings

1. **NEME=1 or 2 is not safe.** DEME=(EFERMI-real(EBOT))/(NEME-1) (aes.f:292–293) for 1 divides by zero. Interpolation clamps index to NEME-1 and uses index-1 and index+1 (rows 513–540). For 2 it can reach the index 0. Therefore, the definition rejects TASK.NEME<3. The same problem occurs if the existing TAU has less than 3 energy, even at native NEME=10.
2. **F core shells may exceed the cache name table.** STRME has only ME1..ME12 (aes.f:79–80), while f shell gives NCST=14. The core state loop uses STRME(ICST) in the absence of a common file ME (row 363). The definition f shells rejects even if an exceptional branch with an existing common ME will bypass the table. This is a conservative limitation of the new definition to the repair of Fortran.
3. **The range of principal shells ends with n=5.** SHELL contains only K,L,M,N,O (aes.f:81). The definition allows n<=5 and controls l<n.
4. **Automatic cache ME is without verification of origin.** Routine itself searches file ME and then MEi, reads matrix elements and does not verify match potential or energy mesh. The existing cache also shuts down PRINTME. The same working directory with previous calculations can change the result even without the input switch. The reader uses cache indexes directly, without a range check. The definition does not manage this native cache.
5. **Capacity TAU is not locally verified.** NE is taken from TAU metadata after allocation according to NEMAX. ETAB, IAES and ETABFIN are then used according to this NE (rows 226–236, 484–486). NE<=NEMAX and the consistency of the input with TAU; the test serialization shall not verify this.
6. **Mesh must be uniform and growing.** DE is taken only from the first two ETAB energy and the whole convolution uses this step. The matrix-element mesh ends at EF while EBOT is the first TAU energy. If EBOT=EF, DEME=0 and for otherwise valid NEME. Source does not verify these assumptions directly in AES.
7. **The output is not a simple table NE rows.** The header is created with NE initial energy but for each core state with writes 2*NE-1 final energy. Each point has two lines and the intensity are complex values. The specialized parser must respect this structure.
8. **The output energy remains in Ry.** aes.f:657–659 writes real(ETABFIN) and the relative difference without RY_EV; the interpretation of eV from spectral output NRAES cannot be transferred here. The final energy shall be calculated from 2*EF-ECOR with the Broadening taken from the first initial energy (lines 483–487).

Numerical accuracy AES, cache behaviour and specific TAU file have not been verified running SPR-KKR. Patch does not perform parser output or does not fix Fortran.

## Tests

13 new cases: registration and defaults, omitted values and EMAX, file roundtrip and pickle; explicit parameters and roundtrip; three hazardous values NEME; three incompatible relativistic regimes; both correlations; invalid shell and limit n; insufficient capacity; f shell exceeding STRME.

Full input_parameter set: **121 passed**, Python 3.12.

```bash
PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 /workspace/scratch/e79ca15d9801/test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q
```

Passed git diff --check and checked the patch reverse application. Numerical SPR-KKR run not started.
