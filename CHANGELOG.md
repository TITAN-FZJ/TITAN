# Changelog

## [Unreleased]


## 2020-09-04

### Added

- Initial directives for openACC

### Fixed

- Fixed a "doubled" loop left behind in the previous changes



## 2020-09-01

### Added

- fullhk is also used for ground state L
- fullhk is also used for time propagation

### Changed

- Removed "hermitianization" of the hamiltonian. In all the tests, there were no differences with and without
- Separated a few do concurrent loops to avoid possible race conditions if the loop is done concurrently (and removed !$OMP directives that were causing errors)
- Removed git tags and increased commit description size from initial git information

### Fixed

- Fixed the calculation of ground state L when supercond=T (it had different results with lambda=0 and supercond=F)



## 2020-08-29

### Changed

- Removed 'addresults' from time-propagation, as 'checkpoints' is a better way to continue the results

### Fixed

- Removed flag '-qpre-fetch=5' that made time-propagation calculations slower
- Fixed local hamiltonian calculation for time-propagation calculation (must be inside time-loop to update the U-term)
- Fixed writing of time-propagated quantities when checkpoint is used



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
