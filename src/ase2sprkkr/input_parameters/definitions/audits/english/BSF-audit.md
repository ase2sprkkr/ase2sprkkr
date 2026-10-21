# BSF — audit of the shared BSF/BSFEK/BSFKK definition

Status: **meaningful with documented defaults and geometry constraints**. This audit describes the current `tasks-workflow` branch and the supplied native SPR-KKR 9.7 source (`PUB9.7_260420_18h13`). It covers the E-K and K-K modes selected through the shared BLOCHSF implementation; it is not a numerical certification of every KPATH/Bravais combination.

## Current ASE2SPRKKR contract

`BSF`, `BSFEK`, and `BSFKK` are public names for one native calculation: `TASK=BSF`, `CONTROL.ADSI=BLOCHSF`. The aliases therefore share one input/output contract rather than having alias-specific filenames.

- `NK1`/`NK2` require at least 2 points. Native BLOCHSF divides by `NK1-1` and `NK2-1` in K-K mode.
- `NK` in E-K mode requires at least 2 points; `NE` requires at least 1.
- `GRID` is mode-dependent: E-K uses 3 by default, K-K uses 0. This avoids the invalid one-point integration situation produced by `GRID=3`, `NE=1`.
- K-K validation rejects non-finite vectors and zero or collinear planes relative to `KA`.
- `K1` and `K2` are endpoints, not direction vectors. Native code subtracts `KA`, then orthogonalizes the second relative direction when needed.

## File contract

The TAU directory used by BLOCHSF is an **optional reusable cache**, not a required input: native code can recompute TAU data when the cache is absent.

Persistent result artifacts are:

- `{dataset}_BLOCHSF.bsf`
- `{dataset}_BLOCHSF_spol.bsf`

The same filenames apply to `BSF`, `BSFEK`, and `BSFKK`, because all three select `ADSI=BLOCHSF`.

Transient scratch files are not part of the task result manifest.

## Energy and grid defaults

### E-K mode

The current explicit defaults are `NE=200`, `GRID=3`, `NK=300`, `ImE=0.001 Ry`, `EMIN=-0.2 Ry`, and `EMAX=-1 Ry`.

`EMAX=-1` is a real absolute energy, **not** a sentinel for the Fermi level. Native `init_mod_energy` substitutes `EF` only for values below approximately `-99`. The default interval therefore decreases from `-0.2` to `-1 Ry`; users should normally set the intended energy window explicitly.

Native reading defaults `NK` to 800, while ASE2SPRKKR explicitly writes 300. This is a deliberate ASE-side default, not a claim of convergence.

### K-K mode

The current explicit energy is `EMIN=0.7 Ry`, with `EMAX=EMIN`, `NE=1`, and `GRID=0`. This is not automatically the Fermi energy. A Fermi-surface calculation must set the intended absolute or relative energy explicitly.

Unset `KA` is completed by native code as the zero origin; ASE2SPRKKR validation uses the same interpretation. `NK1`, `NK2`, `K1`, and `K2` have no invented defaults and must be supplied where required.

## Geometry and native limitations

E-K versus K-K is selected by the actual energy count (`NE>1` versus `NE=1`), not by separate native dispatch keys. ASE2SPRKKR validates conflicts between the requested public alias and the resulting regime.

In E-K mode native code distributes points according to segment lengths, with at least two points per segment. Degenerate explicit segments or zero total path length remain native/geometric edge cases requiring sensible user input.

For 2D systems, the native K-K branch zeros the z components of `K1`/`K2` before subtracting `KA`. A plane that is valid in the supplied 3D vectors can therefore degenerate after native 2D projection. The input definition cannot determine this solely from the vectors; 2D calculations require a geometrically consistent in-plane specification.

BLOCHSF may reuse cached TAU data. ASE2SPRKKR models that reuse but does not prove physical consistency of arbitrary externally supplied cache contents.

Advanced native switches such as `BSFSPN`, `BSFORB`, `BSFETD`, `BSFALF`, and `EPHOT` are outside this audit scope.

## Verification scope

Current regression coverage checks the shared aliases, E-K/K-K regime behavior, K-K `GRID=0`, equality of `EMIN`/`EMAX` in the single-energy regime, round-trip/pickle behavior, point-count minima, finite/nondegenerate geometry, and the shared BLOCHSF result manifest for all three public task names.

No numerical SPR-KKR reference run is claimed by this audit.

Native source references: `blochsf.f` (including the K-K geometry and `NK1-1`/`NK2-1` divisions), `init_mod_energy.f`, and `epath.f`.
