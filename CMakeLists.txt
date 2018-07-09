# Minimum required Version
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

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
                                  source/dtdksub.F90
                                  source/invers.F90
                                  source/U_matrix.F90
                                  source/calc_initial_Uterms.F90
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
                                  source/TCM.F90
                                  source/createFolder.F90
                                  source/eintshechi.F90
                                  source/eintshe.F90
                                  source/calculate_chi.F90
					                        source/calculate_all.F90
				                          source/calculate_dc_limit.F90
                                  source/fermi_surface.F90
                                  source/Atom_variables.F90
                                  source/mod_rotation_matrices.F90
                                  source/check_emergency_stop.F90
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
                                  )

if(PLATFORM MATCHES juqueen)
	set(TITAN_EXE titan_juqueen.exe)
  set(CMAKE_Fortran_COMPILER mpixlf2008_r)
elseif(PLATFORM MATCHES jureca)
	set(TITAN_EXE titan_jureca.exe)
	set(CMAKE_Fortran_COMPILER mpiifort)
elseif(PLATFORM MATCHES osx)
  set(TITAN_EXE titan.exe)
  set(CMAKE_Fortran_COMPILER mpif90)
elseif(PLATFORM MATCHES booster)
  set(TITAN_EXE titan_booster.exe)
  set(CMAKE_Fortran_COMPILER mpiifort)
else()
	set(TITAN_EXE titan.exe)
  set(CMAKE_Fortran_COMPILER mpiifort)
endif()

add_executable(${TITAN_EXE} ${SOURCE_TITAN})

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

if(PLATFORM MATCHES juqueen)
	target_link_libraries(${TITAN_EXE} $ENV{ND}/lib/libnag_essl.a -L/bgsys/local/lib -lesslbg -qsmp=omp)
  link_directories(/bgsys/local/lib)
  target_compile_options(${TITAN_EXE} PRIVATE -WF,-D_JUQUEEN
                                           -O3
                                           -qstrict
                                           -qnoescape
                                           -qnosave
                                           -qsmp=omp
                                           -qthreaded)
  set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
  if(BIGINT)
     target_compile_options(${TITAN_EXE} PRIVATE -qintsize=8)
  endif()
elseif(PLATFORM MATCHES jureca)

  target_link_libraries(${TITAN_EXE} -mkl -lnag -qopenmp)
  target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                           -fpp
                                           -O3
                                           -qopenmp
                                           -mkl
                                           -xHost)
elseif(PLATFORM MATCHES osx)
  target_link_libraries(${TITAN_EXE} -mkl -L$ENV{HOME}/lib -lkibe)
  target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                           -D_OSX
                                           -fpp
                                           -O3
                                           -mkl
                                           -xHost)
elseif(PLATFORM MATCHES booster)
  link_directories(/homeb/jias15/jias1501/BoosterNAG)
  target_link_libraries(${TITAN_EXE} -mkl
                                     -qopenmp
                                     -I/homeb/jias15/jias1501/BoosterNAG/nag_interface_blocks
                                     libnag.a )
  target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                              -fpp
                                              -O3
                                              -mkl
                                              -I/homeb/jias15/jias1501/BoosterNAG/nag_interface_blocks
                                              -xMIC-AVX512
                                              -qopenmp)
elseif(PLATFORM MATCHES "RWTH")
  target_link_libraries(${TITAN_EXE} -mkl $ENV{FLAGS_MATH_LINKER} -qopenmp)
  target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                           -fpp
                                           -O3
                                           -g
                                           -qopenmp
                                           -mkl -$ENV{FLAGS_MATH_INCLUDE}
                                           -xHost)
else()
  target_link_libraries(${TITAN_EXE} -mkl -qopenmp /opt/NAG/fll6i26dcl/lib/libnag_nag.a)
  target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                           -fpp
                                           -O2
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
#					-qopenmp
#					-mkl
#					-xHost
#)

if(PLATFORM MATCHES juqueen)
if(DEBUG)
target_compile_options(${TITAN_EXE} PRIVATE -C
				    	    -g
					    -qfullpath
					    -O0
					    -qflttrap=enable:invalid:overflow:underflow:zerodivide
					    -qinitauto=7FF7FFFF
					    -qnoescape
					    -qkeepparm
					    )
endif(DEBUG)

else()
if(DEBUG)
target_compile_options(${TITAN_EXE} PRIVATE -C -CB
                                        -check uninit

                                        -debug all
                                        -warn all
                                        -ftrapuv
                                        -traceback
                                        -g
                                        -O0
                      )
endif(DEBUG)
endif()
#target_link_libraries(${EXEC_NAME} ${NAG}/lib/libnag_nag.a -qopenmp)
