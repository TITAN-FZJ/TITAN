# Changelog

## [Unreleased]

## 2020-08-28

### Added

- Interface for calcLGS

### Changed

- Changed some variable names of intrinsic Fortran functions

### Fixed

- Variables for f_n and f_n_negative changed to private in calcLGS_eigenstates


## 2020-08-27

### Added

- This CHANGELOG was added (previous changes are listed in commit messages)
- Calculation of Tight-Binding Hamiltonian for all k-points before self-consistency (checks if free memory is enough)
- Added get_memory function in mod_tools to obtain free memory

### Changed

- Supercond hamiltonian is build directly on hamilt_local and hamiltk subroutines
- Green functtion is built directly on green
- Added interfaces and procedure pointers to avoid if(lsupercond) or if(leigenstates)
- Added "only" in many modules (still more to do) to avoid non-desired variable use
- Changed type name System -> System_type to use intrinsic function

### Fixed

- Orbital angular momentum with supercond=T
- Turn off supercond for initial_expectations
- Deallocate local hamiltonian after initial_expectations
- A few missing `.d -> _dp` variables
