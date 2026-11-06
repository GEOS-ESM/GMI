# GMI_GridComp: class modules and the GmiChemistry tree

**Audience:** coding agents and developers working on the GMI chemistry code.
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GMI` @ `ec9e2fd` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Parent documents: [`../docs/ARCHITECTURE.md`](../docs/ARCHITECTURE.md)
(the MAPL component, `GMIchem_GridCompMod.F90`) and the Chemistry document
`../../docs/ARCHITECTURE.md`. This directory
holds no MAPL component of its own; it is the Fortran code that
`GMIchem_GridCompMod` calls. Scientific kernels are described by directory
only.

## 1. Purpose

Run the GMI processes on the species bundles that `GMIchem_GridCompMod`
owns: emissions, forced boundary conditions, deposition, surface area
densities, photolysis, thermal rate constants, chemical solver.

## 2. Place in the hierarchy

`GMI_GridCompMod.F90` defines the derived type `GMI_GridComp` with seven
members, one per class module, and the driver routines
(`GMI_GridCompMod.F90`, type `GMI_GridComp`):

| Member | Class module (file) | Process |
|---|---|---|
| `gcEmiss` | `GmiEmiss_GCCMod` (`GmiEmiss_GridCompClassMod.F90`) | surface, aircraft, lightning, biogenic emissions; also `GmiEmiss_initSurfEmissBundle` |
| `gcFBC` | `GmiForcingBC_GCCMod` | forced surface boundary conditions (reads `<species>_BC` imports) |
| `gcDepos` | `GmiDepos_GCCMod` | dry and wet deposition, gravitational settling |
| `gcSAD` | `GmiSAD_GCCMod` | aerosol surface area densities |
| `gcPhot` | `GmiPhotolysis_GCCMod` | photolysis rates (fast-JX or Cloud-J), aerosol optics |
| `gcThermalRC` | `GmiThermalRC_GCCMod` | thermal rate constants |
| `gcChem` | `GmiChem_GCCMod` (type name `GmiChemistry_GridComp`) | the chemistry solver |

File names end `_GridCompClassMod.F90`; module names end `_GCCMod`. Each
exposes `<Name>_GridComp` (type), `<Name>_GridCompInitialize`, `_GridCompRun`,
`_GridCompFinalize` (public interfaces of each file). They are not
ESMF or MAPL components: no `SetServices`, no own states. They receive the
MAPL import and export states of `GMICHEM` as arguments
(`GMI_GridCompMod.F90`, `GMI_GridCompInitialize`).

## 3. Entry points and call order (`GMI_GridCompMod.F90`)

Common signature (class Run): `(self, bgg, bxx, impChem, expChem, nymd,
nhms, tdt/cdt, ...)`: the two `Species_Bundle` arguments (transported and
not), the import and export `ESMF_State`, date, time and time step.

| Driver | Called from | Order inside |
|---|---|---|
| `GMI_GridCompInitialize` | `Initialize_` | Emiss, Depos, SAD, Photolysis, ForcingBC, ThermalRC, Chemistry Initialize; `GmiEmiss_initSurfEmissBundle`; copies `so4v_saexist` from Photolysis to ThermalRC; reads the `do*` switches |
| `GMI_GridCompRun1` | Run phase 1 | `doEmiss`: Emiss Run (`mixPBL` false); `doBoCo`: ForcingBC Run |
| `GMI_GridCompRun2` | Run phase 2 | `doDepos`: Depos Run (every heartbeat); only when chemistry is due: `doSAD` SAD, `doPhot` Photolysis, `doTherm` ThermalRC (receives `gcSAD` and `gcPhot`), `doSolver` Chemistry |
| `GMI_GridCompRunOrig` | Run phase 99 (single phase) | Depos, Emiss (`mixPBL` true), SAD, Photolysis, ForcingBC, ThermalRC, Chemistry |
| `GMI_GridCompFinalize` | `Finalize_` | Finalize of all seven |

Dependencies between class objects exist only where arguments are passed:
ThermalRC reads SAD and photolysis results, and the shared species bundles
carry results to the solver. *Inferred:* other couplings go through the MAPL
exports (`gmiSAD`, `gmiQJ`, `gmiQK`, ...) and the bundle arrays
(`GmiSwapSpeciesBundlesMod.F90`, `SwapSpeciesBundles`).

## 4. State

- Each class object holds the grid indices `i1 i2 im j1 j2 jm km` and
  lon/lat copies set by `sendToGCs` in `GMIchem_GridCompMod.F90`; deposition
  and photolysis also hold the orbit and clock; the emission class holds
  `heartBeat` and `doingPredictorNow` (`GMIchem_GridCompMod.F90`,
  `sendToGCs`, Run_).
- Each class reads its configuration from `GMI_GridComp.rc` with
  `rcEsmfReadTable` and ESMF config calls in its Initialize routine
  (`GmiShared/GmiESMF`).
- Photolysis state includes `aeroProviderName` and `AM_I_AERO_PROVIDER`
  (set by `Initialize_`) and a fast-JX bundle type from
  `GmiFastJX_includeMod.F90`.
- `GmiSwapSpeciesBundlesMod.F90` provides `SwapSpeciesBundles` and
  `speciesReg_for_CCM` to move species between GEOS bundles and GMI's
  arrays; used by all class modules.
- Class-specific kernels use the `t_Chemistry`, `t_Emission`,
  `t_Deposition`, `t_SpeciesConcentration`, `t_gmiGrid` and `t_GmiClock`
  types from the directories below.

## 5. Configuration that changes structure

- Mechanism: chosen at build time (`GMI_MECHANISM`), see the main document.
- `GMI_GridComp.rc` of the mechanism: `doEmiss`, `doBoCo`, `doDepos`,
  `doSAD`, `doPhot`, `doTherm`, `doSolver`; `phot_opt: 3` with
  `fastj_opt` 4 (fast-JX 6.5, `photolysis/fast_JX65`) or 5 (Cloud-J,
  `photolysis/CloudJ`) (`GmiPhotolysis_GridCompInitialize`; the mechanism
  files set `fastj_opt: 5`); `chem_opt`, `sad_opt`, `do_drydep`,
  `do_wetdep`, `do_simpledep`, `do_emission`, `emissionSpeciesNames::`,
  `extdataBcSpeciesNames::`, `fixedConcentrationSpeciesNames::`
  (`StratTrop_Orig/GMI_GridComp.rc`).
- ExtData: `GMI_ExtData.yaml` names the file collections (monthly emission
  inventories, boundary conditions, SAD, climatologies) and the exports
  (`GMI_GridComp/GMI_ExtData.yaml`, sections `Collections`, `Samplings`,
  `Exports`).

## 6. The `GmiChemistry` tree (directory level)

| Directory | Content | Called by |
|---|---|---|
| top-level `Gmi*_mod.F90` | `GmiChemistryMethod_mod` (type `t_Chemistry`), `GmiUpdateChemistry_mod`, `GmiSolverInterface_mod`, `GmiPhotolysisRateConstants_mod`, `GmiThermalRateConstants_mod`, `GmiUpdateSAD_mod`, `GmiUpdateForcingBC_mod`, `GmiStratosphericLoss_mod`, `GmiShipEmission_mod`, `GmiSurfaceEmissionInChemistry_mod`, `GmiGravitationalSettling_mod`, `GmiSavedVariables_mod`, `GmiFastJX_includeMod` | class modules (`USE` lists) |
| `smv2chem/` | the SMV2 sparse-matrix Gear solver (`smv2chem_solverGMI.F90`, `smvgearGMI.F90`, `jsparseGMI.F90`, `physprocGMI.F90`) | Chemistry class, through the solver interface |
| `photolysis/CloudJ`, `photolysis/fast_JX65` | Cloud-J (7.4) and fast-JX 6.5 photolysis codes; `photolysis/include` shared headers | Photolysis class (`initializeFastJX74`, `initializeFastJX65`) |
| `sad/` | `condense.F90`, aerosol condensation for SAD | SAD class through `GmiUpdateSAD_mod` |
| `AerosolDust/` | `GmiAerDustODSA_mod.F90`, dust optical depth and surface area | *Inferred:* SAD or ThermalRC class |
| `ioChemistry/` | `ReadForcedBC_mod`, `ReadSolarCycle_mod` | ForcingBC and Photolysis classes (`USE` lists) |
| `include/` | `gmi_forc_bc.h`, `smv2chem_par.h` | build include path |
| `<mechanism>/` | per-mechanism generated-style files (below) | build selects one |

Mechanism directories: `Standard`, `Standard_Pyro`, `StratTrop_HFC_S`,
`StratTrop_HFC_S_Pyro`, `StratTrop_JPL19`, `StratTrop_Minimal`,
`StratTrop_Orig` (default). Every one has the same 17 file names:
`GMI_GridComp.rc`, `GMI_Mech_Registry.rc`, `GMI_Registry.rc`,
`setkin_chem_mech.txt`, `setkin_depos.h`, `setkin_fastj.h`,
`setkin_group.h`, `setkin_ibcb.h`, `setkin_kcalc.F90`, `setkin_lchem.h`,
`setkin_mw.h`, `setkin_par.h`, `setkin_ratecalc.F90`, `setkin_smv2.h`,
`setkin_smv2par.h`, `setkin_surf_bc.h`, `setkin_synspc.h` (file lists
compared with `StratTrop_Orig`; only contents differ). Meaning: the
`setkin_*` files are the mechanism definition (species, reactions,
rate calculation); `GMI_Registry.rc` defines the MAPL import, internal and
export specs and feeds the ACG; `GMI_Mech_Registry.rc` lists the species
tables read at run time. *Unverified:* the `setkin_*` files are written by an
external mechanism-setup tool (naming only).

Other directories under `GMI_GridComp/`:

- `GmiEmission/`: `GmiEmissionMethod_mod.F90` (type `t_Emission`),
  `emiss_update.F90`; subdirectories `GOCARTroutines`, `gsfc`, `harvard`,
  `llnl`, `lightning`, `ioEmission` (emission schemes and readers).
- `GmiDeposition/`: `GmiDepositionMethod_mod.F90`, `drydep_update.F90`,
  `wetdep_update.F90`, `simpledep_update.F90`.
- `GmiSpeciesConcentration/spcConcentrationMethod`: species concentration
  type and fixed-concentration handling.
- `GmiShared/`: `GmiESMF` (rc-file reading, field-bundle helpers),
  `GmiIOutilities`, `GmiSupportingModules`, `GmiInclude` (shared headers,
  including `setkin_par.h` consumers' constants).

## 7. Aerosol optics in this layer

Photolysis reads the `AERO` state's `implements_aerosol_optics_method`
attribute and executes `run_aerosol_optics` on it when the aerosol
provider is `GOCART2G` (`GmiPhotolysis_GridCompClassMod.F90`,
`Acquire_OptProps`). The label GMI attaches is described in the
[main document](../docs/ARCHITECTURE.md), section 7.

## 8. Dependencies

`ESMF`, `MAPL`, `Chem_Mod`, `Chem_UtilMod`, `Species_BundleMod` and
`Species_ArrayMod` (`USE` lists of the class modules). Compile-time
mechanism constants via `setkin_par.h` (`NSP`, `NUM_K`, `NUM_J`) read by
`gmi_acg.pl`.

## 9. Known traps

1. Class modules are not MAPL components: do not look for `SetServices`
   there; specs are in `GMIchem_GridCompMod.F90` and the generated headers.
2. The module name and the file name differ (`_GCCMod` vs
   `_GridCompClassMod.F90`); a `grep` for the module name in file names
   finds nothing.
3. Compile-time constants (`NSP`, reaction counts) come from the chosen
   mechanism; array sizes change with the mechanism.
4. `GMI_GridCompRunOrig` and `Run1`/`Run2` use different orders and a
   different `mixPBL` setting for emissions.

## 10. Decisions

Recorded in [`../docs/decisions.md`](../docs/decisions.md).

## 11. Observations (noted, not fixed)

1. `GmiSwapSpeciesBundlesMod.F90` is described in the file header only as
   tools "to transfer species concentrations between GEOS-5 and GMI"; the
   detailed copy semantics were not read.
2. Open question for maintainers: `GmiChemistry/AerosolDust` and `sad/`
   callers were inferred from names, not from `USE` lists.
