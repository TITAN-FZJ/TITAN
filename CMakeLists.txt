# Minimum required Version
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

# Projects name
project(TITAN)

option(DEBUG "Enable debugging" OFF)

enable_language(Fortran)

file(GLOB sources source/*.F90
	  	  source/*.f90)
set(SOURCE_TITAN ${SOURCE_TITAN}  source/main.F90
                                  source/create_files.F90
                                  source/check_files.F90
                                  source/sort_all_files.F90
                                  source/ldos.F90
                                  source/ldos_and_coupling.F90
                                  source/ldos_energy.F90
                                  source/ldos_jij_energy.F90
                                  source/calculate_chi.F90
                                  source/fermi_surface.F90
                                  source/Npl_variables.F90
                                  source/green.F90
                                  source/hamiltk.F90
                                  source/invers.F90
                                  source/U_matrix.F90
                                  source/sumkshechi.F90
                                  source/build_rotation_matrices.F90
                                  source/eintshechi.F90
                                  source/L_gs.F90
                                  source/check_emergency_stop.F90
                                  source/band_structure.F90
                                  source/mod_f90_kind.F90
                                  source/mod_constants.F90
                                  source/mod_parameters.F90
                                  source/mod_mpi_pars.F90
                                  source/mod_tools.F90
                                  source/mod_input.F90
                                  source/mod_dnsqe.F90
                                  source/mod_progress.F90
                                  source/EnergyIntegration.F90
                                  source/AtomTypes.F90
                                  source/mod_System.F90
                                  source/mod_monoBasis.F90
                                  source/mod_polyBasis.F90
                                  source/Lattice.F90
                                  source/mod_BrillouinZone.F90
                                  source/mod_SOC.F90
                                  source/SK_TightBinding.F90
                                  source/TightBinding.F90
                                  source/mod_magnet.F90
                                  source/ElectricField.F90
                                  source/mod_io.F90
                                  source/mod_self_consistency.F90
                                  source/mod_susceptibilities.F90
                                  source/mod_alpha.F90
																)
add_executable(titan.exe ${SOURCE_TITAN})

if(PLATFORM MATCHES "Juqueen")
  set(CMAKE_Fortran_COMPILER mpixlf2008_r)
else()
  set(CMAKE_Fortran_COMPILER mpiifort)
endif()

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

if(PLATFORM MATCHES "Juqueen")
  link_directories(/bgsys/local/lib)
  target_link_libraries(titan.exe $ENV{ND}/lib/libnag_essl.a -lesslbg -qsmp=omp)
  target_compile_options(titan.exe PRIVATE -WF,-D_JUQUEEN
                                           -O3
                                           -qstrict
                                           -qnoescape
                                           -qnosave
                                           -qsmp=omp
                                           -qthreaded)
  set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
elseif(PLATFORM MATCHES "Jureca")

  target_link_libraries(titan.exe -mkl -lnag -qopenmp)
  target_compile_options(titan.exe PRIVATE -no-wrap-margin
                                           -fpp
                                           -O3
                                           -qopenmp
                                           -mkl
                                           -xHost)
else()
  target_link_libraries(titan.exe -mkl -qopenmp /opt/NAG/fll6i26dcl/lib/libnag_nag.a)
  target_compile_options(titan.exe PRIVATE -no-wrap-margin
                                           -fpp
                                           -O3
                                           -qopenmp
                                           -mkl
                                           -xHost)
endif()









#target_compile_options(${EXEC_NAME} PRIVATE -no-wrap-margin
#                                        -fpp
#                                        -D _LINUX
#                                        -O3
##                                        -I${NAG}/c_headers#
#					-qopenmp
#					-mkl
#					-xHost
#)
if(DEBUG)
target_compile_options(titan.exe PRIVATE -CB
                                        -check uninit
                                        #-check-all
                                        -debug all
                                        -warn all
                                        -ftrapuv
                                        -traceback
                                        -g
                                        -O0
                      )
endif(DEBUG)

#target_link_libraries(${EXEC_NAME} ${NAG}/lib/libnag_nag.a -qopenmp)
