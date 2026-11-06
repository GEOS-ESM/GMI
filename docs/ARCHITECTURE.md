# GMI architecture

**Audience:** coding agents and developers working on the GMI repository.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GMI` @ `fce93fe` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

**Build status on this line:** not built, not ported. On `release/MAPL-v3`, `@GEOSchem_GridComp/CMakeLists.txt` calls `esma_add_subdirectories(GOCART)` only, with the comment "Restore subdirs as components are ported to MAPL3."; `GMI` appears only in a commented-out list of subdirectories. The Chemistry top-level library (`GEOS_ChemGridComp.F90`, which adds this component) is present but not compiled on this line: the `esma_add_library` block that builds it is commented out in `../CMakeLists.txt`. Statements about it describe the file as it is. This repository's `CMakeLists.txt` lists the dependencies `Chem_Shared` (from the Chemistry `Shared` directory, also in the commented-out list), `GEOS_Shared`, `esmf`, `Process_Library`. The code here is the MAPL2 code of an older state of the repository, kept for a later port; it is one commit ("Updates for MAPL3", which renames `use MAPL_ConstantsMod` to `use MAPL_Constants` in six files under `photolysis/CloudJ` and `GmiEmission`, and strips trailing whitespace) ahead of, and 20 commits behind, MAPL2 `main`.

**Differences from the MAPL2 line.** This line lacks later MAPL2 changes: the `Standard_Pyro` mechanism directory; the later `Standard` mechanism update (species `CH3CN`, `HCN`, ODS emission fields); the `FLOOR_VALUE`, `Br_limit`, `HBr_limit`, `BrONO2_limit` run-time keys and the Br spike traps; the debug exports `CPU_MAP`, `SO2_MMR`, `SO2CMASS`; and tendencies for non-transported species. Sections below describe the code on this line.

GMI is the Global Modeling Initiative combined stratosphere/troposphere
chemistry package, wrapped as one MAPL gridded component, `GMICHEM`. The
GEOS interface layer is `GMIchem_GridCompMod.F90`; the scientific code is
under `GMI_GridComp/`. This document covers the interface layer in depth
and the scientific code at directory level. The class-module layer and the
`GmiChemistry` tree are described in
[`GMI_GridComp/ARCHITECTURE.md`](../GMI_GridComp/ARCHITECTURE.md).

## 1. Placement

- Fixture path: `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@GMI`
  (nested repository, `components.yaml` key `GMI`).
- Parent: the Chemistry component,
  `../../docs/ARCHITECTURE.md` (it adds `GMICHEM`, runs it and builds the `AERO`
  state). It adds the child with `MAPL_AddChild(GC, NAME='GMICHEM',
  SS=GMI_SetServices)` when `ENABLE_GMICHEM:` is true; `GMI_SetServices` is
  the `use ... only: GMI_SetServices => SetServices` rename of
  `GMIchem_GridCompMod` (`GEOS_ChemGridComp.F90`, SetServices).
- Children: none in the MAPL sense. `GMICHEM` has no `MAPL_AddChild`; its
  seven sub-processes are plain Fortran derived types called directly
  (`GMI_GridComp/GMI_GridCompMod.F90`, type `GMI_GridComp`).
- Build target: one library defined by `CMakeLists.txt`; not built on this line. The parent links it
  as `GMI` (the parent's build).

## 2. Directory and file map

| Path | What it is |
|---|---|
| `GMIchem_GridCompMod.F90` | The MAPL gridded component: `SetServices`, Initialize, Run1/Run2, Finalize, `AERO` state, aerosol optics method. Only public symbol: `SetServices`. |
| `CMakeLists.txt`, `gmi_acg.pl` | Build; selects one chemical mechanism; runs the GMI-specific code generator before MAPL's ACG. |
| `GMI_GridComp/GMI_GridCompMod.F90` | Aggregate type `GMI_GridComp` holding the seven class objects; `GMI_GridCompInitialize`, `Run1`, `Run2`, `RunOrig`, `Finalize`. |
| `GMI_GridComp/Gmi{Chem,Depos,Emiss,ForcingBC,Photolysis,SAD,ThermalRC}_GridCompClassMod.F90` | The class modules (module names `Gmi*_GCCMod`): one process each. |
| `GMI_GridComp/GmiSwapSpeciesBundlesMod.F90` | Copies species between the MAPL species bundles and GMI's arrays. |
| `GMI_GridComp/GMI_ExtData.rc`, `GMI_ExtData.yaml` | ExtData (MAPL external data) specification for GMI's file inputs; installed to `etc`. |
| `GMI_GridComp/GmiChemistry/` | Chemistry solver, photolysis, SAD code and the per-mechanism directories (`Standard`, `StratTrop_HFC_S`, `StratTrop_HFC_S_Pyro`, `StratTrop_JPL19`, `StratTrop_Minimal`, `StratTrop_Orig`). |
| `GMI_GridComp/GmiEmission/`, `GmiDeposition/`, `GmiSpeciesConcentration/`, `GmiShared/` | Emission, deposition, species-concentration and shared support code. |
| `README.md`, `CHANGELOG.md`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `LICENSE*`, `ChangeLog` | Repository metadata. `README.md` is a one-line title. |

There was no `docs/` directory before these documents.

## 3. Purpose (component `GMICHEM`)

`GMICHEM` advances the species of one GMI chemical mechanism: emissions,
forced boundary conditions, dry and wet deposition, surface area densities,
photolysis, thermal rate constants and the chemical solver. It exports
oxidants and species to GOCART, ACHEM, MAM, TR and CARMA, and optionally
provides the `AERO` state (prescribed aerosols, with optics) to Radiation
(`GMIchem_GridCompMod.F90`; `GEOS_ChemGridComp.F90`, SetServices).

## 4. Place in the hierarchy

- Parent: Chemistry (`GEOS_ChemGridComp`), child name `GMICHEM`.
- Parent-side wiring (all in `GEOS_ChemGridComp.F90`, SetServices):
  - `CHEMENV` to `GMICHEM`: `AIRDENS`, `DELP`, `LFR`, `TPREC`,
    `LIGHT_NO_PROD`.
  - `GMICHEM` to `ACHEM`: `OH`, `H2O2`, `NO3`, `OX` (as `O3`); also `O`,
    `OH`, `OX`, `OCS_JRATE` as `O3P`, `OHSTRAT`, `O3`, `OCS_JRATE` when
    both are enabled.
  - `GMICHEM` to `GOCART` (`OH`, `CH4`, `H2O2`, `NO3`, `HNO3`; `OX` as
    `O3`), to `GOCART2G` (as `GMI_OH` and so on), to `MAM` (`ALK4`, `CH4`,
    `C2H6`, `C3H6`, `C3H8`, `ISOP`), to `TR` (`OX` as `OX_TR`, `stOX_loss`,
    `DD_OX`) and to `CARMA` (`H2SO4`, `HNO3`).
  - `AERO` from `GOCART2G` or `GOCART.data` into `GMICHEM` according to the
    provider name; `CARMA_SUSAREA`, `CARMA_SUREFF` to `SO4SAREA`, `SO4REFF`
    when the provider is `CARMA`; `HEMCO` `GMI_ISOPRENE` when
    `doMEGANviaHEMCO` is true.
  - Parent sets `SULFURIC_ACID_SOURCE:` to `full_field` when CARMA and GMI
    are both on; GMI reads it in Run.
- Friendly transport: transported species are `FriendlyTo`
  `DYNAMICS:TURBULENCE:MOIST` plus the component name, so they land in the
  parent's `CHEM_TRACERS` bundle. `OX` is also friendly to `ANALYSIS` when
  `ANALYSIS_OX_PROVIDER:` is `GMICHEM` (`GMIchem_GridCompMod.F90`,
  SetServices).

## 5. Entry points (`GMIchem_GridCompMod.F90`)

- **SetServices**: creates two `Runtime_Registry` objects from
  `GMI_Mech_Registry.rc` (tables `GMI_table::` transported, `XX_table::`
  not transported); registers Initialize, `Run1`, `Run2` (two Run phases)
  and Finalize with `MAPL_GridCompSetEntryPoint`; stores a private
  `GMIchem_state` as internal state `GMIchem_state`; adds import, internal
  and export specs; calls `MAPL_GenericSetServices`.
- **Initialize_**: `MAPL_GenericInitialize`; `extract_` (clock, time step,
  state); sets photolysis's `aeroProviderName` and `AM_I_AERO_PROVIDER`;
  `sendToGCs` copies grid bounds, orbit, clock and lon/lat into the seven
  class objects; creates the two `Species_Bundle` objects with no
  allocation and points them at the MAPL internal-state arrays; calls
  `GMI_GridCompInitialize`; sets up overpass masks (10 am and 2 pm local
  time, `OVP_init`); builds the `AERO` state if GMI is the provider;
  initialises chemical families. On this line it reads no species floor or
  Br trap keys from `GMI_GridComp.rc`.
- **Run1** (phase 1): if `GMI_RUN_PHASES:` is 2 (default) calls `Run_` with
  phase 1; otherwise does nothing.
- **Run2** (phase 2): calls `Run_` with phase 2, or phase 99 when
  `GMI_RUN_PHASES:` is 1.
- **Run_**: reads the run alarm (`RunGMINow`, chemistry time step set by
  `GMICHEM_DT:`, default the heartbeat; error if shorter); applies a fixed floor of 1.0e-30 to the transported species (`#define FLOOR_VALUE`); stores
  initial species to form tendencies (transported species only); loads `OCS_CLIMO`; sets `H2SO4` to 1e-30
  unless `SULFURIC_ACID_SOURCE:` begins `full_field`; then
  - phase 1: `GMI_GridCompRun1` (emissions then forced boundary
    conditions, every heartbeat);
  - phase 2: `GMI_GridCompRun2` (deposition every heartbeat; SAD,
    photolysis, thermal rate constants, solver only when `RunGMINow`),
    then turns the run alarm off;
  - phase 99: `GMI_GridCompRunOrig` (all processes in the original order,
    only when `RunGMINow`);
  then fills the tendency exports (for example `H2O_TEND`) and the overpass
  fields (`OVP10_*`, `OVP14_*`) (`GMIchem_GridCompMod.F90`, Run1, Run2,
  Run_). *Unverified:* the rest of `Run_` after phase 99 was only skimmed.
- **Finalize_**: `GMI_GridCompFinalize`, destroys bundles, registries and
  masks, `MAPL_GenericFinalize`.

The phase numbers visible to the parent are the registration order; the
parent's phase handling is in the parent's architecture document.

## 6. State

- **Imports**: declared in `SetServices` by hand (`AERO` per provider,
  `SO4SAREA`/`SO4REFF` for CARMA, `<species>_BC` from the
  `extdataBcSpeciesNames::` table in `GMI_GridComp.rc`, `SHIP_NO` if
  `do_ShipEmission`, `RI`/`RL` if `fastj_opt` is 5, `OCS_CLIMO`, `CNV_FRC`,
  `GMI_ISOPRENE`, `LIGHT_NO_PROD`, `VEG_FRAC`, `LAI_FRAC`) plus the
  generated include `GMICHEM_ImportSpec___.h` from the mechanism's
  `GMI_Registry.rc`.
- **Internal state**: one spec per registry species (transported and not)
  plus the families `Bry` and `Cly` (restart skipped). Non-transported species
  have `ADD2EXPORT` (`GMIchem_GridCompMod.F90`, SetServices). Species
  restart is the MAPL internal restart of `GMICHEM`; the parent template
  keeps the `GMICHEM` restart-file keys commented out (see the parent's
  architecture document). `importRestartFile: gmichem_import_rst` is set in each
  mechanism's `GMI_GridComp.rc`.
- **Exports**: hand-written (`AERO`,
  `AERO_ACI`, `AERO_DP`, `HNO3CONDsad`, `HNO3GASsad`, `jNO2val`, bundles
  `surfEmissForChem`, `gmiSAD`, `gmiQJ`, `gmiQQJ`, `gmiQK`, `gmiQQK`,
  `gmiERADIUS`, `gmiTAREA`, overpass `OVP10_*`/`OVP14_*`, `stOX_loss`) plus
  generated includes `GMICHEM_ExportSpec___.h`, `Deposition_ExportSpec___.h`,
  `Reactions_ExportSpec___.h`, `Tendency_ExportSpec___.h`. `AERO`,
  `AERO_ACI`, `AERO_DP` exist only when `AERO_PROVIDER:` is `GMICHEM` or
  `CARMA`.
- **Species ownership**: GMI species live in MAPL's internal state;
  `Species_Bundle` objects (`bgg`, `bxx`) hold pointers into it
  (`Initialize_`, `MAPL_GetPointer` on the internal state). The class
  modules obtain working copies through `SwapSpeciesBundles`
  (`GMI_GridComp/GmiSwapSpeciesBundlesMod.F90`).

## 7. The `AERO` state and its methods

- Built in `Aero_StateInitialize` (contained in `Initialize_`) only when
  `AERO_PROVIDER:` is `GMICHEM` (`AM_I_AERO_PROVIDER`). It fills `AERO` with
  a bundle `AEROSOLS` of 15 prescribed aerosols (`BCphobic`, `BCphilic`,
  `du001`-`du005`, `OCphobic`, `OCphilic`, `ss001`-`ss005`, `SO4`; these are
  GMI exports declared in `GMI_Registry.rc` as "prescribed ... from
  GMICHEM"; the photolysis class fills them from the same-named imports
  divided by air density, `GmiPhotolysis_GridCompClassMod.F90`, Run), a Mie table instance (2, created with
  `Chem_MieCreate`), the attribute names for pressure, relative humidity,
  extinction, single-scattering albedo and asymmetry parameter, and sets
  `implements_aerosol_optics_method` true.
- Method label attached: `run_aerosol_optics` (`ESMF_MethodAdd`, user routine
  `aerosol_optics`). It converts mass mixing ratios to column mass per layer
  and calls `mie_` with the GOCART Mie table, for one band.
- `AERO_ACI` gets `implements_aerosol_activation_properties_method` false.
- Photolysis executes the same label on an `AERO` state when the provider is
  `GOCART2G` (`GmiPhotolysis_GridCompClassMod.F90`, `Acquire_OptProps`:
  reads `implements_aerosol_optics_method`, sets `band_for_aerosol_optics`
  per wavelength, `ESMF_MethodExecute(label='run_aerosol_optics')`).
  Cross-repository label story: the parent's architecture document (`AERO`).
- With `GMICHEM` as the provider and `GOCART2G` aerosols, photolysis
  stops with an error (`GmiPhotolysis_GridCompClassMod.F90`, consistency
  check in Initialize).

## 8. Configuration that changes structure

| Key | File | Effect |
|---|---|---|
| `GMI_MECHANISM` | CMake cache option, default `StratTrop_Orig` | Selects the mechanism directory; must have `setkin_chem_mech.txt` (`CMakeLists.txt`). |
| `ENABLE_GMICHEM:` | Chemistry configuration (`GEOS_ChemGridComp.F90`, SetServices) | Adds the child. |
| `AERO_PROVIDER:` | `AGCM.rc` (component config) | `GOCART.data`, `GOCART2G`, `GMICHEM`, `CARMA`, `none`; selects `AERO` import or export; anything else is an error (`SetServices`). |
| `ANALYSIS_OX_PROVIDER:` | `AGCM.rc` | `GMICHEM` makes `OX` friendly to `ANALYSIS`. |
| `GMI_RUN_PHASES:` | `AGCM.rc` | 1 or 2 run phases (default 2). |
| `GMICHEM_DT:` | `AGCM.rc` | Chemistry time step (default heartbeat). |
| `SULFURIC_ACID_SOURCE:` | `AGCM.rc` | `full_field` keeps GMI `H2SO4`; else reset each step. |
| `doEmiss`, `doBoCo`, `doDepos`, `doSAD`, `doPhot`, `doTherm`, `doSolver` | `GMI_GridComp.rc` | Switch individual sub-processes (default true) (`GMI_GridCompInitialize`). |
| `phot_opt`, `fastj_opt` (4 or 5) | `GMI_GridComp.rc` | Photolysis scheme: 4 calls `initializeFastJX65`, 5 calls `InitializeFastJX74` from `CloudJ` (`GmiPhotolysis_GridCompInitialize`). |
| `doMEGANviaHEMCO:`, `do_ShipEmission:` | `GMI_GridComp.rc` | Adds imports (see State). |
| `REPLAY_MODE:` | `AGCM.rc` | Replay predictor step passed to emission class. |

## 9. Build structure

- `CMakeLists.txt` globs sources from a fixed list of directories, plus the
  chosen mechanism directory; the library depends on `Chem_Shared`,
  `GEOS_Shared`, `esmf`, `Process_Library`.
- The mechanism is chosen at CMake configure time, not at run time: files
  of the same names exist in every mechanism directory (the file lists are
  identical) and only the chosen directory is compiled, and its
  `GMI_Registry.rc`, `GMI_Mech_Registry.rc` and `GMI_GridComp.rc` are
  installed to `etc`.
- Code generation: a custom command runs `gmi_acg.pl -R` (extra registry
  files for deposition, reactions, tendencies, from `GMI_Registry.rc`,
  `setkin_lchem.h`, `setkin_depos.h`, `setkin_par.h`), then MAPL's
  `mapl_acg.pl` on each registry (producing `*_ExportSpec___.h`,
  `*_GetPointer___.h`, `*_DeclarePointer___.h`), then `gmi_acg.pl -v`
  (second-stage `GetPointer2` files), then `mapl_acg.pl` on `GMI_Registry.rc`
  with name `GMICHEM` (import/export includes and `GMICHEM_History___.rc`).
  The generated headers are copied into the build's `GMI_GridComp`
  directory.
- Tests: none in this repository (no test directory found).

## 10. Dependencies

Uses `ESMF`, `MAPL`, `Chem_Mod`, `Chem_UtilMod`, `Chem_GroupMod`,
`Species_BundleMod`, `Runtime_RegistryMod`, `OVP` (overpass masks),
`Chem_Mie` from the Chemistry shared code in the parent repository
(`GMIchem_GridCompMod.F90`, `USE` list). `CMakeLists.txt` lists `Chem_Shared`, `GEOS_Shared`, `esmf`, `Process_Library`; it does not list `Chem_Base`.

## 11. Known traps

1. **Mechanism is a build choice.** The species list, import/export
   includes and rate arrays change with `GMI_MECHANISM`; a run directory
   needs the `GMI_Mech_Registry.rc` and `GMI_GridComp.rc` of the same
   mechanism.
2. **`GMI_GridComp.rc` is read in several places** (`SetServices`,
   `Initialize_`, `GMI_GridCompInitialize`, each class Initialize, and by
   the parent when HEMCO is on), so a key must be consistent for all
   readers.
3. **Phase numbers.** Phase 99 means "single phase" and is internal to GMI;
   in single-phase mode `Run1` does nothing and `Run2` does everything.
4. **Tendencies.** Phase 1 and 99 store the tendency, phase 2 adds to it
   (comment in `Run_`). On this line tendency exports are restricted to transported species (`Run_` comment; `gmi_acg.pl` keeps a tendency export only when the friends field is `D,T,C`).
5. **Generated headers** (`*___.h`) are not in the repository; searching
   the source tree will not find them.

## 12. Decisions

See [`decisions.md`](decisions.md).

## 13. Observations (noted, not fixed)

1. `SetServices` adds the `AERO`/`AERO_ACI`/`AERO_DP` exports when
   `AERO_PROVIDER:` is `CARMA`, but `Initialize_` fills them only when it is
   `GMICHEM`; a commented-out variant of the condition references CARMA
   (`GMIchem_GridCompMod.F90`, SetServices, Initialize_). Question for
   maintainers: is the `CARMA` export intentionally left empty?
2. In `SetServices` the `AERO_PROVIDER` case `"GMICHEM"` adds no import (its
   specs are commented out); the same value is accepted without an `AERO`
   import.
3. `CMakeLists.txt` lists photolysis source directories (`fast_JX`,
   `fastj`, `fast_JX53b`, `fast_JX53c`, `fastJX53c_ref`, `lookup`, `utils`)
   that are absent from `GMI_GridComp/GmiChemistry/photolysis/` (only
   `CloudJ`, `fast_JX65`, `include` exist), and refers to
   `GMI_GridComp/GmiChem_GridCompMod.F90`, which does not exist (the file
   is `GmiChem_GridCompClassMod.F90`). Question for maintainers: are these
   leftovers?
4. `GMI_GridComp.rc` in `StratTrop_Orig`, `StratTrop_Minimal` and other
   mechanisms writes `doEmiss    T` and `doBoCo     T` without a colon,
   while the code reads the labels `doEmiss:` and `doBoCo:`
   (`GMI_GridCompInitialize`); *Inferred:* the lines are therefore not read
   and the defaults (true) apply.
5. The script header of `gmi_acg.pl` defaults to a directory
   `GmiChemistry/stratTropMech` that does not exist; the CMake call always
   passes explicit paths.
6. `README.md` is one line; no other documentation exists in the
   repository.
