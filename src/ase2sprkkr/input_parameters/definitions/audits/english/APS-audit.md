# APS audit and new definition

Sources: supplied SPR-KKR 9.7 and currently fetched ASE2SPRKKR refactor, commit 267512c4b8ccbbfd37af7ede898febc066b7e3ee. APS does not have its own definition. Single patch 0014-aps.patch adds aps.py and test_aps.py against this base; previous tasks does not include.

## Input and defaults

APS uses kkrgen and NRAPS. initvar.f:281–287 sets INITELOOP=false, RUNELOOP=false, USELPLS1=true and RDDOS=true. kkrmain.f:264 compiles DOSFIL=DATASET.dos and opens this file at RDDOS. NRAPS expects old spin-resolved format DOS with specific header structure. The name of any DOS file does not read this branch from TASK.

| Parameter |Native behaviour|Definitions| Serialization |
| --- | --- | --- | --- |
| TASK.IT | 1 | 1 | default omitted |
| TASK.CL | 2P | 2P | default omitted |
| TASK.CALCLS | false | false | omitted |
| TASK.OUTGNU | false | false | omitted |
| TASK.GAMMA |Locally -0,1; absence shuts off Broadening| None | omitted |
| MODE |IREL=min(1,IREL), at normal original IREL=3 i.e. 1| SREL |explicitly entered|
| ENERGY.GRID/NE/ImE |task does not choose own mesh, energy is later read from DOS| 1 / 300 / 0 Ry |explicitly entered|
| ENERGY.EMIN/EMAX |EMIN for APS is set to EF; then the table is rewritten from DOS|not in this definition| omitted |

GAMMA=None expresses absence and off Broadening, not physical width -0,1. Omitted defaults IT, CL and switches correspond to the source. SREL explicitly expresses the result native switch, NREL remains allowed. At FREL or SP-SREL the definition warns that the source will change the mode.

GRID=1, NE=300 and ImE=0 are explicit definitions, not claims of native APS defaults. NE determines the allocated capacity NEMAX; the actual amount of energy and ETAB are taken from DOS. Better than fixed default 300 would be to take the capacity from the header of a specific file, including, where appropriate, the breaking reserve. This patch automatic reading DOS metadata does not implement.

GAMMA is written as eV and transferred by division RY_EV. Its actual use is not constant Lorentz width: in NRAPS multiplys absolute real energy, see below. CALCLS triggers additional (l,s)-distinguished spectral outputs, OUTGNU additional gnuplot data.

## Findings in the Source

1. **Unprotected capacity DOS.** NRAPS first allocates by NEMAX from ENERGY.NE. Later it reads NE from DOS header and writes ETAB(IE,1) without control NE<=NEMAX (nraps.f:835–865). I without breaking must ENERGY.NE sufficient for the whole file. IT, NT and NL must also correspond to the potential and structure DOS; only input definition does not verify their conformity.
2. **At least two energies and an evenly growing mesh.** Routine without condition reads ETAB(2)-ETAB(1) and continues to work with a single step DE. DOS with one energy is not applicable; non-uniform or descending mesh is not adequately verified. The definition rejects ENERGY.NE<2, but the actual number and mesh of the file must verify the reader DOS.
3. **Broadening may exceed multiple fields.** For GAMMA>0 is added 5 eV below EF. NADD=NINT((5/RY_EV)/DELE)+1; the source controls only NADD+1<=101. Later, it sets NE=NE+NADD and expands ETAB, which has only NEMAX rows (nraps.f:947–954). In addition, NEBRD=2*NADD+NE+1 must enter NEMAX+101. Checking alone NADD is not sufficient. Settings NE must include a reserve; selecting a larger capacity does not replace control of the source.
4. **The width is related to absolute energy.** For E>=EFIN is used WLORTAB=GAMMAV*DREAL(EBRD), not GAMMAV*(E-EFIN) (nraps.f:929–931). Energies below EFIN receive zero. Dependence of result at absolute energy start and unit coefficient deserves physical correction or explanation; parameter cannot be honestly described as normal constant width in eV.
5. **The direction is switched to hard.** initvar.f:779 performed IREL=min(1,IREL); NRAPS stopped in addition IREL>1 (nraps.f:270). FULLPOT zakazuje check_input.f:65–68. DMFT and LDA+U prohibits GEN. The definition rejects both correlations; FULLPOT must capture SPR-KKR from potential.
6. **MPI branch APS does not complete.**APS is not included in the list tasks continuing after MPI energy loop in gen.f:1355–1368. The definition has mpi=False.

Definition validates physical shell l<n and positive 1 IT. The upper limit IT depends on the potential. GAMMA shall be final; hodnota <=0 shall indicate that Broadening is switched off. No GAMMA restrictions derived from Broadening are introduced.

## Tests and verification limits

8 new cases: Discovery and defaults, omitted TASK values, file roundtrip and pickle, explicit shell/type/switches/breading and their roundtrip, incorrect shell, both unsupported correlations, mode switching alerts and off breaking, insufficient energy capacity.

Full input_parameter set: **116 passed**, Python 3.12, 5,36 s.

```bash
PYTHONPATH=src ASE2SPRKKR_NO_USER_PROFILE=1 /workspace/scratch/e79ca15d9801/test-env/bin/python -m pytest src/ase2sprkkr/input_parameters/test -q
```

Passed git diff --check and checked the patch reverse application. Numerical APS run or check a specific DOS file have not been performed. Patch does not contain specialized parser APS outputs.
