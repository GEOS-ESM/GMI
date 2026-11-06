# AGENTS.md

GMI: Global Modeling Initiative combined stratosphere/troposphere
chemistry, one MAPL gridded component (`GMICHEM`) with Fortran class
modules and several selectable chemical mechanisms.

**Placement:** GEOSgcm fixture, `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@GMI`
(`components.yaml` @ `0f2e0fb`, key `GMI`). Parent: Chemistry, see
[`../docs/ARCHITECTURE.md`](../docs/ARCHITECTURE.md) (adds `GMICHEM`, wires
`AERO`). Verified against `GMI` @ `ec9e2fd`; describes `develop`.

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
- Built as part of GEOSgcm through the Chemistry parent; configure with
  `-DGMI_MECHANISM=<name>` to change mechanism (one of the directories under
  `GMI_GridComp/GmiChemistry` holding `setkin_chem_mech.txt`).
- No test suite in this repository.
- Run-time files from the chosen mechanism (`GMI_GridComp.rc`,
  `GMI_Registry.rc`, `GMI_Mech_Registry.rc`) and `GMI_ExtData.*` install to
  `etc`.
