# AIPES separate definitions and audits

AIPES has a dedicated definition `aipes.py` on `tasks-workflow` and uses the `kkrspec` executable. Native-physics assessment: **works but has major caveats** © native AIPES has a specific local/uncoherent approximation meaning that cannot be confused with the numerical integration of all ARPES angles. The Numerical Solver hasn't been launched.

## Actual native branch

SPEC_RSLAB recognizes AIPES and sets IBLOCH=2. In SPEC_TAU_PREPARE for IBLOCH=2 is enforced CPAATOM=1 for all Layer atoms, even if the material is not a chemical alloy. This is the internal branch of the calculation of projectors and incoherent members, not the change of physical occupancy.

SPEC_MECALC calls ZMATINCO for all these atoms. SPEC_SUMUP_INT for AIPES zeros ROTRLA/ROTELA and consists ROALL of ROAT+ROSUR plus ROINC. In contrast, for ARPES, it uses ROAT+ROTRA+ROTER+ROSUR and other correction of non-coherent members. Native therefore explicitly denotes angularly integrated photoemis with surface contribution and implements an approximation branch other than the simple sum of the angular network.

NT/NP are not the number of knots of any general hemisphere quadrature. Here is audited procedure with fixed final electron geometry and NT=NP=1. The definition rejects angular/momentum scans as outside its audited range; it is not a claim that native all their variants are prohibited.

## Definitions and validation

It shares structural and spectroscopic sections ARPES but has its own TASK AIPES, CONTROL ADSI, validator and upper limit of energy. The Photon section remains: unlike SPLEED, they have EPHOT and polarization physical significance.

Fixed final electron angles, NT/NP, absence of unaudited scans moments, non-negative final workfunction and control, positive final EPHOT and one positive energy mesh without SPLITSS are validated. INPVER=0 uses external old input.inp and modern parameters are not considered authoritarian in this case; content and existence of external files this input-only validator does not verify.

## Defaults and take-over difference ARPES

EMINEV=−8 eV and new EMAXEV=0 eV** face EF*** are written to the input. For valence photoemisi the range to EF is more suitable than taken +5 eV, which includes vacant states. The change is explicit ASE option, not unwritten default, and the user can change it. SPEC_INPUT photoemission branches generate their own energy scan; this range cannot be mechanically compared with defaults general GEN contours.

Listed NE=300, EWORK_EV=4.2 eV, IMV_INI_EV=0.05 eV, IMV_FIN_EV=2 eV, EPHOT=25 eV, NL=4, N_LAYER=50, NLAT_G_VEC=57 and N_LAYDBL={10,10}, practical options inherited from ARPES remain. Native photoemission start contains EMESH=200, VIL/VIH=0.1/1 eV, N_LAYER=30 and NLAT_G_VEC=30. These are listed deviations: they cannot be derived from a mapping error. A larger number of layers can be better for convergence, while generic workfunction can be worse than the value of a particular potential.

Electron THETA/PHI=0..0 and NT/NP=1 correspond to photoemission default. SPOL=2 remains explicit photoemission default, versus SPLEED native SPOL=1. Expert polarization conventions have not been fully audited; we do not fake equivalence of missing keys and writing out the same numerical value.

CRYS_VECS=False is ** not written** and leaves the use of primitive crystal vectors. Native presence CRYS_VECS sets USE_CRY_PRIM_VECS=false. The original name CRYS_VEC is not read, therefore the separate definition offers a corrected name and description; the actual `aipes.py` definition must be used for its task-specific contract.

## Next ALE

- Local/uncoherent approximation and surface member should be taken into account when interpreting and comparing with experimental angle-integrated PES. This definition does not demonstrate that each geometry or material gives the same result as the integral ARPES.
- It is necessary to converge layers, reciprocal cutoffs, angular moments, energy mesh and control. Output work and photon energy set according to material and experiment.
- STRVER/INPVER and the surface geometry generated shall correspond to the potential. Full-potential, multiparticle and all spin-density-matrix branches have not been certified throughout the chain.
- mpi=False new definition is conservative choice, not evidence broken native MPI implementation.
- Native double typo LANZ1 assignment and absolute IMV_INI to VIH applies here. Listed N_LAYDBL sets both directional values and used IMV_INI_EV writes correctly VIL.
