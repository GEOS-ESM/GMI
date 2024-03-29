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
