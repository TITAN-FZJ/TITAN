# Changelog

## [Unreleased]

======================================================================================================================

## 2020-09-18

### Changed

- Moved 'build_hext' to 'mod_hamiltonian'
- Moved calls to 'vector_potential' and 'magnetic_field' outside (to t loop)
- Moved calculation of 'hamilt_nof = ( eval_kn(n,ik) * id ) - hk' outside (to mod_time_propagator)
- 'td_hamilt' and the respective procedures were not necessary anymore
- Added if( lfullhk ) to initial call to hamiltk
- Changed KronProd from subroutine to function
- Simplified LS_Solver

======================================================================================================================

## 2020-09-15

### Added

- Added flag '-align array64byte' to try to improve vectorisation

### Changed

- Changed some divisions by respective multiplications

### Fixed

- Time-propagated energy is now being reduced by MPI too
- Time-propagated OAM and singlet_coupling are now being reduced by MPI too

======================================================================================================================

## 2020-09-15

### Added

- Declarations of external functions to avoid warnings on new ifort version

### Changed

- Removed remaining 'do concurrents' as its behaviour still seems to be discussed
- Moved 'invers' subroutine to mod_tools
- Changed 'green.F90' to a module 'mod_greenfunction.F90'
- Changed 'hh' to 'hlrs' in CMakeLists

### Fixed

- Expectation value of the energy in time-propagation was not using the pointer

======================================================================================================================

## 2020-09-11

### Changed

- Removed some do concurrents due to complaints of some compilers

======================================================================================================================

## 2020-09-07

### Fixed

- Fixed time propagation when printfieldonly is used (files were being created after)

======================================================================================================================

## 2020-09-04

### Added

- Initial directives for openACC

### Fixed

- Fixed a "doubled" loop left behind in the previous changes

======================================================================================================================

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

======================================================================================================================

## 2020-08-29

### Changed

- Removed 'addresults' from time-propagation, as 'checkpoints' is a better way to continue the results

### Fixed

- Removed flag '-qpre-fetch=5' that made time-propagation calculations slower
- Fixed local hamiltonian calculation for time-propagation calculation (must be inside time-loop to update the U-term)
- Fixed writing of time-propagated quantities when checkpoint is used

======================================================================================================================

## 2020-08-28

### Added

- Interface for calcLGS

### Changed

- Changed some variable names of intrinsic Fortran functions

### Fixed

- Variables for f_n and f_n_negative changed to private in calcLGS_eigenstates

======================================================================================================================

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

======================================================================================================================