# AGENTS.md

GMI: Global Modeling Initiative combined stratosphere/troposphere
chemistry, one MAPL gridded component (`GMICHEM`) with Fortran class
modules and several selectable chemical mechanisms.

**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GMI` @ `fce93fe` (2026-11-06).
**Placement:** GEOSgcm fixture, `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@GMI`
(`components.yaml` @ `f50622d`, key `GMI`). Parent: Chemistry, see
`../docs/ARCHITECTURE.md` (adds `GMICHEM`, wires
`AERO`).

**Build status on this line:** not built, not ported. On `release/MAPL-v3`, `@GEOSchem_GridComp/CMakeLists.txt` calls `esma_add_subdirectories(GOCART)` only, with the comment "Restore subdirs as components are ported to MAPL3."; `GMI` appears only in a commented-out list of subdirectories. The Chemistry top-level library (`GEOS_ChemGridComp.F90`, which adds this component) is present but not compiled on this line: the `esma_add_library` block that builds it is commented out in `../CMakeLists.txt`. Statements about it describe the file as it is. This repository's `CMakeLists.txt` lists the dependencies `Chem_Shared` (from the Chemistry `Shared` directory, also in the commented-out list), `GEOS_Shared`, `esmf`, `Process_Library`. The code here is the MAPL2 code of an older state of the repository, kept for a later port; it is one commit ("Updates for MAPL3", which renames `use MAPL_ConstantsMod` to `use MAPL_Constants` in six files under `photolysis/CloudJ` and `GmiEmission`, and strips trailing whitespace) ahead of, and 20 commits behind, MAPL2 `main`.

**Differences from the MAPL2 line.** This line lacks later MAPL2 changes: the `Standard_Pyro` mechanism directory; the later `Standard` mechanism update (species `CH3CN`, `HCN`, ODS emission fields; the `Standard` mechanism files differ in content); the `FLOOR_VALUE`, `Br_limit`, `HBr_limit`, `BrONO2_limit` run-time keys and the Br spike traps; the debug exports `CPU_MAP`, `SO2_MMR`, `SO2CMASS`; and tendencies for non-transported species. Sections below describe the code on this line.

## Documents
- [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md): component entry points, state, `AERO`, configuration, build.
- [`GMI_GridComp/ARCHITECTURE.md`](GMI_GridComp/ARCHITECTURE.md): the seven class modules and the `GmiChemistry` tree.
- [`docs/decisions.md`](docs/decisions.md): decisions recorded so far.
- `README.md`, `CHANGELOG.md`, `CONTRIBUTING.md`: repository metadata (not edited by the architecture documents).

## Map
- `GMIchem_GridCompMod.F90`: the MAPL component (`SetServices`, `Run1`, `Run2`, `AERO`).
- `GMI_GridComp/GMI_GridCompMod.F90`: calls the class modules `Gmi*_GridCompClassMod.F90`.
- `GMI_GridComp/GmiChemistry/<mechanism>/`: mechanism files; one is compiled (`GMI_MECHANISM`, default `StratTrop_Orig`).
- `gmi_acg.pl`, `CMakeLists.txt`: generate import/export specs from the mechanism registry.

## Rules for agents
- Note bugs and oddities in the "Observations" sections; do not fix them
  unless asked.
- Keep the documents current when you change entry points, resource
  keys, mechanism choice or generated files.
- Cite `(file, routine)` or `(file, KEY)`; no line numbers. Mark
  unchecked statements *Inferred:* or *Unverified:*.
- Do not read the numerical kernels (`smv2chem/`, `photolysis/`) unless the
  task needs them; the documents stop at directory level there.
- Generated headers (`GMICHEM_ImportSpec___.h`, `*_ExportSpec___.h`) are not in the
  repository; they appear in the build directory.

## Build and test
- Not built on this line (see build status above). Meant to be built as part of GEOSgcm through the Chemistry parent; configure with
  `-DGMI_MECHANISM=<name>` to change mechanism (one of the directories under
  `GMI_GridComp/GmiChemistry` holding `setkin_chem_mech.txt`).
- No test suite in this repository.
- Run-time files from the chosen mechanism (`GMI_GridComp.rc`,
  `GMI_Registry.rc`, `GMI_Mech_Registry.rc`) and `GMI_ExtData.*` install to
  `etc`.
