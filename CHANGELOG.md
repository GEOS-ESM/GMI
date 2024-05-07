# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Fixed
### Added
### Changed
### Removed
### Deprecated


## [1.3.0] - 2024-05-07

### Fixed

- Set the reaction rates for certain ICE and NAT and STS reactions to zero when HCl < min-concentation; this helps prevent HCl problems.

### Added

- Capability to use 2D ExtData files as Boundary Conditions
- Added surface flux for CH2BR2 and CHBR3 in the HFC+S mechanism
- Added an automatic Br adjustment, based on the species present in the mechanism. CH3Br and CH2Br2 are changed, either by adding to Forced BC's or by scaling Emissions.

### Changed

- Updated CI to use new v2 orb
- Changed the Monthly Wrap sampling (in ExtData yaml) to hold values constant for the month (zero diff, since DMS currently does not use that sampling)
- Modified the starting date for several collections in ExtData yaml; now they cover a wider time range
- Extended SO2 fire flux (ExtData yaml) to include future years from SSP2-4.5
- Changed Forced Boundary Conditions to use the ExtData implementation by default (ASCII still supported); RefD1 and RefD2 BCs available via ExtData, use RefD2 by default.
- Now automatically scale CH2Br2 emissions by 1.8, when using the HFC+S mechanism, to account for missing Br
- Changed the tolerance for the SMV Gear solver; this helps with convergence


## [1.2.0] - 2023-12-07

### Fixed

- QQJK diag calc had bug, used time(t+dt) constituents ibstead of time(t) fields, fixed on 2023Jul19
- Fixed bug in Emissions and Deposition where SZA degrees was used instead of cosine(SZA)
- Fixed units for degassing volcano point emissions
- Improved handling of phot_opt cases; allows for phot_opt == 0 again
- Improved handling of sad_opt cases; allows for sad_opt == 0 again

### Added

- Capability for AERO_PROVIDER=GOCART2G data-driven (GOCART.data)
- Export OCS_JRATE for use by ACHEM

### Changed

- Update CI to use Baselibs default from the CircleCI orb
- Changed a few SAD reactions to use JPL 2019 approach
- ExtData now uses Benchmark G configuration
- Emissions including 2D, 3D and point source (volcano) now are done as in Benchmark G
- Isoprene scaling changed from 1 to 0.7
- Boundary condition file now includes HFCs, and no longer includes extra Br (when HFC mech is chosen)


## [1.1.0] - 2023-04-24

### Added

- Capability for AERO_PROVIDER=GOCART2G
- Capability for AERO_PROVIDER=GMICHEM
- Capability for AERO_PROVIDER=CARMA (uses GMICHEM right now)

### Changed

- Instead of calling the GMI routines for SZA (in Chem_Shared), now call a wrapper for the corresponding MAPL routine
- Changed default aerosdust filename (for case AERO_PROVIDER=GMICHEM)

### Removed

- Code that might have eventually run GOCART-like aerosols from within the GMI framework (note: doubtful this ever worked and would take significant effort to get working IMHO)

## [1.0.0] - 2023-01-18

### Added

- Initial release of GMI based on GMIchem_GridComp code from GEOSchem_GridComp v1.11.0
