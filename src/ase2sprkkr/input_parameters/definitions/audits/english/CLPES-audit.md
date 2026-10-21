# CLPES — current definition audit

`clpes.py` exposes the native KKRSPEC core-level photoemission branch with a conservative serial execution policy. Native restrictions below remain relevant; no numerical accuracy claim is made.
## Parameters and defaults
SPEC_CL is the actual name of the section (spec_inputs.f:1877–1921). Implementation IT (explicitly 1), CL (explicitly 2P), IE_ICST (default 1). CL is always written: without the present CL reader leaves 1S, although SECTION_SET_STRING states fallback 2P. Optional CE is not exposed in this definition; IT determines potential type and native selection includes also same atomic count Z (spec_clpes.f:213–216), not isolated site selection.

NE=300 written, single mesh without SPLITSS. Driver rewrites the energy window in the range of core energy extended by +/-3 eV (spec_clpes.f:301–306); EMIN/EMAX are not user spectral limits. NE=1 selected IE_ICST and validated 4*l+2. Core energy is related to EF (r.224).

EPHOT=1000 eV is explicitly writes instead of inherited 25 eV: more suitable for core PES, but the user must verify the threshold of a specific shell. This is not about guaranteeing the availability of each core level. Divided explicit WORK=4.2, IMV_INI_EV=.05, IMV_FIN_EV=2, NL=4, N_LAYER=50, NLAT_G_VEC=57 are retained; are not declared native defaults. CRYS_VEC corrected for actual CRYS_VECS; its presence turns off primitive crystal vectors. NPOL absent corresponds to native 1; 0 and 3 initiates STOP and rejects them.

## Major caveats
- NCPA>0 vede k STOP (spec_clpes.f:157). Definition does not verify occupation of potential and cannot guarantee the ordered system.
- The audited range is a fixed final electron geometry. Angular and moment scans are rejected; old external input INPVER=0 their control goes around.
- NE=1 uses ECOR/XNORM by loop of atoms (spec_clpes.f:310–319); for multiple types/site may correspond to the last processed set, not the desired type. This branch cannot be considered to be generally verified.
- The selected shell must exist as core shell, IT must exist and EPHOT must exceed the threshold. Python validation without potential will not guarantee these conditions.
- MPI has not been certified, automatic use is disabled; this is not a claim that native MPI is missing.

Category for future joint audit: works but has a large ALE. It's not a setholder off.
