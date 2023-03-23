# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed

### Fixed

### Removed

### Deprecated


### Added

- Capability for AERO_PROVIDER=GMICHEM
- Capability for AERO_PROVIDER=CARMA (uses GMICHEM right now)
- Connectivity from GOCART2G to GMI chem

### Removed

- Code that might have eventually run GOCART-like aerosols from within the GMI framework (note: doubtful this ever worked and would take significant effort to get working IMHO)


## [1.0.0] - 2023-01-18

### Added

- Initial release of GMI based on GMIchem_GridComp code from GEOSchem_GridComp v1.11.0
