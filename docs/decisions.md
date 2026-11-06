# GMI: previous architectural decisions

**Audience:** coding agents and developers working on the GMI repository.
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GMI` @ `ec9e2fd` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

### Chemical mechanism chosen at build time
- Pattern: `GMI_MECHANISM` selects one directory under
  `GMI_GridComp/GmiChemistry`; only it is compiled, and its registry and rc
  files are installed (`CMakeLists.txt`, `GMI_MECHANISM`, `ALLOWED_MECHANISM`). [Verified]
- Likely reason: species count and reaction count are compile-time
  constants (`setkin_par.h`) and the MAPL specs are generated from them. [Inferred]
- Question for maintainers: is run-time mechanism selection ruled out by
  the array sizes only?

### Import and export specs are generated
- Pattern: `gmi_acg.pl` and MAPL's `mapl_acg.pl` generate
  `GMICHEM_ImportSpec___.h`, `GMICHEM_ExportSpec___.h` and the deposition,
  reaction and tendency variants from the mechanism registry
  (`CMakeLists.txt`, `first_stage_files`; `GMIchem_GridCompMod.F90`,
  `#include` lines in SetServices). [Verified]
- Likely reason: the species and reaction lists differ by mechanism. [Inferred]

### Processes are plain Fortran classes inside one MAPL component
- Pattern: seven class objects are called in a fixed order by driver
  routines; none is a MAPL child (`GMI_GridComp/GMI_GridCompMod.F90`). [Verified]
- Likely reason: the code predates the MAPL child structure (module
  revision history starts 2003). [Inferred]

### Two run phases, with a one-phase fallback
- Pattern: phase 1 does emissions and boundary conditions, phase 2
  deposition and (at the chemistry time step) the rest; `GMI_RUN_PHASES:` 1
  collapses them into phase 99 in the original order
  (`GMIchem_GridCompMod.F90`, Run1, Run2, Run_; `GMI_GridCompMod.F90`,
  `GMI_GridCompRunOrig`). [Verified]
- Likely reason: let the parent split emission from chemistry around other
  components, keeping the old behavior available. [Inferred]
- Question for maintainers: which order is used in production runs?

### Chemistry time step separate from the heartbeat
- Pattern: `GMICHEM_DT:` (default the heartbeat) with a run alarm; deposition
  runs every heartbeat while SAD, photolysis, rate constants and the
  solver run when the alarm rings (`GMIchem_GridCompMod.F90`, `extract_`,
  Run_; `GMI_GridCompRun2`). [Verified]

### Aerosols prescribed, optics through a method label
- Pattern: when `AERO_PROVIDER:` is `GMICHEM`, GMI builds an `AERO` state
  of 15 prescribed aerosols and attaches `run_aerosol_optics`
  (`GMIchem_GridCompMod.F90`, `Aero_StateInitialize`). [Verified]
- Likely reason: GMI aerosols are data-driven (comment in the routine: "GMI
  aerosols are prescribed"). [Verified comment]
