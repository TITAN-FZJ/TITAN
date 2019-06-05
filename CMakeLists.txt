# Minimum required Version
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Projects name
project(TITAN)

option(DEBUG "Enable debugging" OFF)
option(BIGINT "Use 8byte integer" OFF)

enable_language(Fortran)

set(SOURCE_TITAN ${SOURCE_TITAN}  source/main.F90
                                  source/AtomTypes.F90
                                  source/Lattice.F90
                                  source/create_files.F90
                                  source/check_files.F90
                                  source/sort_all_files.F90
                                  source/green.F90
                                  source/hamiltk.F90
                                  source/mod_distributions.F90
                                  source/dtdksub.F90
                                  source/invers.F90
                                  source/mod_U_matrix.F90
                                  source/mod_initial_expectation.F90
                                  source/mod_expectation.F90
                                  source/initConversionMatrices.F90
                                  source/build_U.F90
                                  source/mod_System.F90
                                  source/mod_polyBasis.F90
                                  source/mod_BrillouinZone.F90
                                  source/mod_SOC.F90
                                  source/mod_magnet.F90
                                  source/mod_LDOS.F90
                                  source/mod_Coupling.F90
                                  source/ElectricField.F90
                                  source/mod_io.F90
                                  source/ldos.F90
                                  source/ldos_and_coupling.F90
                                  source/ldos_energy.F90
                                  source/ldos_jij_energy.F90
                                  source/coupling.F90
                                  source/jij_energy.F90
                                  source/mod_TCM.F90
                                  source/createFolder.F90
                                  source/eintshechi.F90
                                  source/eintshe.F90
                                  source/calculate_chi.F90
                                  source/calculate_all.F90
                                  source/calculate_dc_limit.F90
                                  source/fermi_surface.F90
                                  source/mod_Atom_variables.F90
                                  source/mod_rotation_matrices.F90
                                  source/mod_check_stop.F90
                                  source/band_structure.F90
                                  source/mod_f90_kind.F90
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
                                  source/SK_TightBinding.F90
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
                                  source/read_band_points.F90
                                  source/TorqueTorqueResponse.F90
                                  source/TorqueSpinResponse.F90
                                  source/mod_RK_matrices.F90
                                  source/mod_imRK4.F90
                                  source/mod_time_propagator.F90
                                  source/mod_imRK4_parameters.F90
                                  )

# Execulable names
if(PLATFORM MATCHES jureca)
  set(TITAN_EXE titan_jureca.exe)
  set(CMAKE_Fortran_COMPILER mpiifort)
elseif(PLATFORM MATCHES osx)
  set(TITAN_EXE titan_osx.exe)
  set(CMAKE_Fortran_COMPILER mpif90)
elseif(PLATFORM MATCHES booster)
  set(TITAN_EXE titan_booster.exe)
  set(CMAKE_Fortran_COMPILER mpiifort)
elseif(PLATFORM MATCHES iff)
  set(TITAN_EXE titan_iff.exe)
  set(CMAKE_Fortran_COMPILER mpiifort)
elseif(PLATFORM MATCHES hh)
  set(TITAN_EXE titan_hh.exe)
  set(CMAKE_Fortran_COMPILER ftn)
else()
  set(TITAN_EXE titan.exe)
  set(CMAKE_Fortran_COMPILER mpiifort)
endif()

add_executable(${TITAN_EXE} ${SOURCE_TITAN})

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

execute_process(COMMAND git describe --abbrev=4 --dirty --always --tags OUTPUT_VARIABLE GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

if(PLATFORM MATCHES jureca)
  target_link_libraries(${TITAN_EXE} -mkl -lnag -qopenmp)
  target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                              -fpp
                                              -DVERSION="${GIT_VERSION}"
                                              -qopenmp
                                              -mkl
                                              -xHost)
elseif(PLATFORM MATCHES osx)
  target_link_libraries(${TITAN_EXE} -mkl -L$ENV{HOME}/lib -lkibe)
  target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                              -D_OSX
                                              -fpp
                                              -DVERSION="${GIT_VERSION}"
                                              -qopenmp
                                              -mkl
                                              -xHost)
elseif(PLATFORM MATCHES booster)
  link_directories($ENV{HOME}/BoosterNAG)
  target_link_libraries(${TITAN_EXE} -mkl
                                     -qopenmp
                                     -I$ENV{HOME}/BoosterNAG/nag_interface_blocks
                                     -L$ENV{HOME}/BoosterNAG
                                     -lnag )
  target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                              -fpp
                                              -DVERSION="${GIT_VERSION}"
                                              -mkl
                                              -I$ENV{HOME}/BoosterNAG/nag_interface_blocks
                                              -xMIC-AVX512
                                              -qopenmp)
elseif(PLATFORM MATCHES RWTH)
  target_link_libraries(${TITAN_EXE} -mkl $ENV{FLAGS_MATH_LINKER} -qopenmp)
  target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                              -fpp
                                              -DVERSION="${GIT_VERSION}"
                                              -g
                                              -qopenmp
                                              -mkl -$ENV{FLAGS_MATH_INCLUDE}
                                              -xHost)
elseif(PLATFORM MATCHES iff)
  target_link_libraries(${TITAN_EXE} -mkl -L$ENV{HOME}/lib -lkibe)
  target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                              -D_OSX
                                              -DVERSION="${GIT_VERSION}"
                                              -fpp
                                              -static-intel
                                              -mkl
                                              -xHost)
elseif(PLATFORM MATCHES hh)
  target_link_libraries(${TITAN_EXE} -mkl -L$ENV{HOME}/lib -lnag -qopenmp)
  target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                              -fpp
                                              -D_OSX
                                              -DVERSION="${GIT_VERSION}"
                                              -qopenmp
                                              -mkl
                                              -xCORE-AVX2)
else()
  target_link_libraries(${TITAN_EXE} -mkl -qopenmp /opt/NAG/fll6i26dcl/lib/libnag_nag.a)
  target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                           -fpp
                                           -DVERSION="${GIT_VERSION}"
                                           -qopenmp
                                           -mkl
                                           #-ipo
                                           #-qopt-report=5
                                           -g
                                           -xHost)
endif()

#target_compile_options(${EXEC_NAME} PRIVATE -no-wrap-margin
#                                        -fpp
#                                        -D _LINUX
#                                        -O3
##                                        -I${NAG}/c_headers#
#         -qopenmp
#         -mkl
#         -xHost
#)

if(DEBUG)
  target_compile_options(${TITAN_EXE} PRIVATE -C
                                              -CB
                                              -check uninit
                                              -debug all
                                              -warn all
                                              -ftrapuv
                                              -traceback
                                              -qopenmp
                                              -fpp
                                              -DVERSION="${GIT_VERSION}"
                                              -g
                                              -O0)
else()
  target_compile_options(${TITAN_EXE} PRIVATE -O3)
endif(DEBUG)
#target_link_libraries(${EXEC_NAME} ${NAG}/lib/libnag_nag.a -qopenmp)
