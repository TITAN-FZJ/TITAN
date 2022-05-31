# Minimum required Version
cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
# Generates a 'compile_commands.json' file containing the exact compiler calls for all translation units of the project in machine-readable form.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Since GNU Compiler is always available, look first for mpiifort which is usually linked with intel
find_program(CMAKE_Fortran_COMPILER NAMES $ENV{FC} mpiifort mpif90 PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{FC} mpicxx PATHS ENV PATH NO_DEFAULT_PATH)

# Projects name
project(TITAN Fortran CXX)

# Defining options
option(DEBUG "Enable debugging" OFF)
option(PREP "Enable Instrumentation" OFF)

# Source files
set(SOURCE_TITAN ${SOURCE_TITAN}  source/main.F90
                                  source/AtomTypes.F90
                                  source/Lattice.F90
                                  source/create_files.F90
                                  source/check_files.F90
                                  source/sort_all_files.F90
                                  source/mod_greenfunction.F90
                                  source/mod_hamiltonian.F90
                                  source/mod_distributions.F90
                                  source/dtdksub.F90
                                  source/mod_U_matrix.F90
                                  source/mod_initial_expectation.F90
                                  source/mod_expectation.F90
                                  source/mod_System.F90
                                  source/mod_polyBasis.F90
                                  source/mod_BrillouinZone.F90
                                  source/mod_SOC.F90
                                  source/mod_magnet.F90
                                  source/mod_LDOS.F90
                                  source/mod_coupling.F90
                                  source/ElectricField.F90
                                  source/mod_io.F90
                                  source/mod_logging.F90
                                  source/coupling.F90
                                  source/jij_energy.F90
                                  source/mod_TCM.F90
                                  source/createFolder.F90
                                  source/eintshechi.F90
                                  source/eintshe.F90
                                  source/calculate_chi.F90
                                  source/calculate_all.F90
                                  source/calculate_dc_limit.F90
                                  source/mod_fermi_surface.F90
                                  source/mod_Atom_variables.F90
                                  source/mod_rotation_matrices.F90
                                  source/mod_check_stop.F90
                                  source/mod_band_structure.F90
                                  source/mod_kind.F90
                                  source/mod_constants.F90
                                  source/mod_parameters.F90
                                  source/mod_mpi_pars.F90
                                  source/mod_tools.F90
                                  source/mod_input.F90
                                  source/mod_dnsqe.F90
                                  source/mod_chkder.F90
                                  source/mod_progress.F90
                                  source/mod_sumrule.F90
                                  source/EnergyIntegration.F90
                                  source/TightBinding.F90
                                  source/mod_self_consistency.F90
                                  source/mod_susceptibilities.F90
                                  source/mod_alpha.F90
                                  source/mod_beff.F90
                                  source/mod_torques.F90
                                  source/mod_disturbances.F90
                                  source/mod_prefactors.F90
                                  source/adaptiveMesh.F90
                                  source/setLoops.F90
                                  source/TorqueTorqueResponse.F90
                                  source/TorqueSpinResponse.F90
                                  source/mod_superconductivity.F90
                                  source/mod_time_propagator.F90
                                  source/mod_cuda.F90
                                  source/mod_nvtx.F90
                                  source/mod_dft.F90
                                  )

# Setting instrumentation of the code when -DPREP=on
set(INSTRUMENT)
if(PREP)
  set(CMAKE_Fortran_COMPILER scorep-mpif90)
  set(CMAKE_CXX_COMPILER scorep-nvcc)
  set(INSTRUMENT _inst)
  message("Instrumentation on")
  message("Don't forget to export the line: export SCOREP_WRAPPER_INSTRUMENTER_FLAGS='--thread=omp'")
endif()

# Getting the COMPILER when not given
if(NOT COMPILER)
  set(COMPILER ${CMAKE_Fortran_COMPILER_ID})
  message("Setting COMPILER=${COMPILER}")
endif()
# Changing $COMPILER variable to lowercase for comparison later
string(TOLOWER ${COMPILER} COMPILER)


# Checking preparation/instrumentalization flag
if(NOT PREP)
  set(PREP OFF)
endif()

# Adding debugging suffix when present
if(DEBUG)
  set(DEBUGSUFFIX _debug)
  message("DEBUG on")
else()
  set(DEBUGSUFFIX)
endif()

# Adding compiler suffix:
# INTEL
if(COMPILER MATCHES intel)
  set(CMAKE_Fortran_COMPILER mpiifort)
  set(COMPILERSUFFIX _intel)
  message("Using Intel Compiler")
endif()
# GNU
if(COMPILER MATCHES gfortran OR COMPILER MATCHES gcc OR COMPILER MATCHES gfort OR COMPILER MATCHES gnu)
  set(COMPILER gnu)
  set(CMAKE_Fortran_COMPILER mpif90)
  set(COMPILERSUFFIX _gnu)
  message("Using GNU Compiler")
endif()
# NVIDIA
if(COMPILER MATCHES pgi OR COMPILER MATCHES nv OR COMPILER MATCHES nvidia)
  set(COMPILER nv)
  set(COMPILERSUFFIX _nv)
  message("Using NVIDIA Compiler")
  set(CMAKE_CXX_COMPILER nvcc)
  set(SOURCE_TITAN_CPP ${SOURCE_TITAN_CPP}  source/diagonalize_gpu.cpp
                                  )
endif()

# Setting coverage suffix
set(COVSUFFIX)
if(COV AND COMPILER MATCHES gnu)
  set(COVSUFFIX _cov)
endif(COV AND COMPILER MATCHES gnu)

# Adding other suffixes when present
if(SUFFIX)
  set(SUFFIX _${SUFFIX})
else()
  set(SUFFIX)
endif()

# Defining executable name
set(TITAN_EXE titan${DEBUGSUFFIX}${COMPILERSUFFIX}${INSTRUMENT}${COVSUFFIX}${SUFFIX}.exe)

# Add an executable to the project using the specified source files
add_executable(${TITAN_EXE} ${SOURCE_TITAN} ${SOURCE_TITAN_CPP})
# Set (main) language
set_property(TARGET ${TITAN_EXE} PROPERTY LINKER_LANGUAGE Fortran)
# Set location for binary
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Add libraries
find_package(LAPACK REQUIRED)
find_package(OpenMP)
target_link_libraries(${TITAN_EXE} ${LAPACK_LIBRARIES} OpenMP::OpenMP_Fortran)

# Getting GIT version
execute_process(COMMAND git describe --abbrev=8 --dirty --always OUTPUT_VARIABLE GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

if(PREP)
  #target_link_libraries(${TITAN_EXE} --thread=omp) 
  target_compile_options(${TITAN_EXE} PRIVATE -D_OLDMPI ) 
endif(PREP)
#target_compile_options(${TITAN_EXE} PRIVATE -D_OLDMPI ) 

# Setting coverage flags (only for GNU compiler)
if(COV)
  if(COMPILER MATCHES gnu)
    message("Coverage on")
    target_link_libraries(${TITAN_EXE} -coverage -fprofile-abs-path) 
    target_compile_options(${TITAN_EXE} PRIVATE -coverage -fprofile-abs-path) 
  else(COMPILER MATCHES gnu)
    message("WARNING: Coverage not activated! Please use -DCOMPILER=gnu")
  endif(COMPILER MATCHES gnu)
endif(COV)

# Setting DEBUG definition
if(DEBUG)
  target_compile_definitions(${TITAN_EXE} PRIVATE DEBUG)
endif(DEBUG)

if(COMPILER MATCHES intel) # INTEL:
  if(NOT DEBUG)
    target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                                -fpp
                                                -DVERSION="${GIT_VERSION}"
                                                -g
                                                -qopt-zmm-usage=high
                                                -align array64byte
                                                -march=core-avx2
                                                -O3
                                                )
  else(NOT DEBUG)
    target_compile_options(${TITAN_EXE} PRIVATE -C
                                                -CB
                                                -check uninit
                                                -debug all
                                                -warn all
                                                -e18
                                                -ftrapuv
                                                -traceback
                                                -fpp
                                                -DVERSION="${GIT_VERSION}"
                                                -g
                                                -O0
                                                )
  endif(NOT DEBUG)
elseif(COMPILER MATCHES gnu) # GFORTRAN:
  if(NOT DEBUG)
    target_compile_options(${TITAN_EXE} PRIVATE -cpp
                                                -DVERSION="${GIT_VERSION}"
                                                -fallow-argument-mismatch
                                                -ffree-form
                                                -ffree-line-length-0
                                                -g
                                                -O3
                                                -march=native
                                                )

  else(NOT DEBUG)
    target_compile_options(${TITAN_EXE} PRIVATE -C
                                                -Wall
                                                -Wextra
                                                -Wconversion
                                                -Wline-truncation
                                                -Wcharacter-truncation
                                                -Wsurprising
                                                -Waliasing
                                                -Wunused-parameter
                                                -Wno-maybe-uninitialized
                                                -fwhole-file
                                                -fimplicit-none
                                                -fallow-argument-mismatch
                                                -fbacktrace
                                                -ffree-form
                                                -ffree-line-length-0
                                                -fcheck=all
                                                -ffpe-trap=zero,overflow,underflow
                                                -finit-real=nan
                                                -ftrapv
                                                -cpp
                                                -DVERSION="${GIT_VERSION}"
                                                -g
                                                -O0
                                                )
  endif(NOT DEBUG)

elseif(COMPILER MATCHES nv) # NVIDIA:
  target_compile_definitions(${TITAN_EXE} PRIVATE _GPU )
  # target_link_libraries(${TITAN_EXE} -L$ENV{MKLROOT}/lib/intel64 -Wl,-rpath,$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -I${CUDA_ROOT}/include -Mcuda -L${CUDA_ROOT}/lib64 -lcusolver -lcudart -lnvToolsExt -acc -gpu=cc70,cc80)

  set_target_properties(${TITAN_EXE} PROPERTIES LINK_FLAGS "-I${CUDA_ROOT}/include -Mcuda -L${CUDA_ROOT}/lib64 -lcusolver -lcudart -lnvToolsExt -acc -gpu=cc70,cc80")
  if(NOT DEBUG)
    string(CONCAT FORTRAN_FLAGS "-cpp "
                                "-DVERSION=\\\"${GIT_VERSION}\\\"  "
                                "-lineinfo "
                                "-O3 "
                                "-lm "
                                "-fast "
                                "-acc "
                                "-gpu=cc70,cc80 "
                                "-Minfo=accel "
                                "-Mcuda "
                                "-g "
                                )
    set_source_files_properties(${SOURCE_TITAN} PROPERTIES COMPILE_FLAGS ${FORTRAN_FLAGS} )
    string(CONCAT CXX_FLAGS "--gpu-architecture=sm_80 "
                            "-lineinfo "
                            "-g "
                            )
    set_source_files_properties(${SOURCE_TITAN_CPP} PROPERTIES COMPILE_FLAGS ${CXX_FLAGS} )

  else(NOT DEBUG)
    string(CONCAT FORTRAN_FLAGS "-cpp "
                                "-DVERSION=\\\"${GIT_VERSION}\\\" "
                                "-lineinfo "
                                "-g "
                                "-Mcuda "
                                "-Mdclchk "                          # instructs the compiler to require that all program variables be declared
                                "-Minform=inform "                   # instructs the compiler to display all error messages (inform, warn, severe and fatal).
                                "-Mchkptr "                          # instructs the compiler to check for pointers that are dereferenced while initialized to NULL (Fortran only).
                                "-Mchkstk "                          # instructs the compiler to check the stack for available space in the prologue of a function and before the start of a parallel region.
                                "-traceback "
                                "-Ktrap=fp,unf "                     # fp, unf is equivalent to inv,divz,ovf, unf  correspond to the processor’s exception mask bits: invalid operation, divide-by-zero, overflow, and underflow, respectively.
                                "-gpu=cc70,cc80 "
                                "-Minfo=accel "
                                "-acc "
                                "-O0 "
                                )
    set_source_files_properties(${SOURCE_TITAN} PROPERTIES COMPILE_FLAGS ${FORTRAN_FLAGS})
    string(CONCAT CXX_FLAGS "--gpu-architecture=sm_80 "
                            "-g "
                            "-G "
                            )
    set_source_files_properties(${SOURCE_TITAN_CPP} PROPERTIES COMPILE_FLAGS ${CXX_FLAGS} )
  endif(NOT DEBUG)
endif(COMPILER MATCHES intel)

# unset(DEBUG CACHE)
# unset(DEBUGSUFFIX CACHE)
# unset(PREP CACHE)
# unset(SUFFIX CACHE)
# unset(COMPILER CACHE)
