# Minimum required Version
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Projects name
project(TITAN Fortran CXX)

option(DEBUG "Enable debugging" OFF)
option(PREP "Enable Instrumentation" OFF)
option(BIGINT "Use 8byte integer" OFF)

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
                                  source/initConversionMatrices.F90
                                  source/mod_System.F90
                                  source/mod_polyBasis.F90
                                  source/mod_BrillouinZone.F90
                                  source/mod_SOC.F90
                                  source/mod_magnet.F90
                                  source/mod_LDOS.F90
                                  source/mod_coupling.F90
                                  source/ElectricField.F90
                                  source/mod_io.F90
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
                                  source/mod_RK_matrices.F90
                                  source/mod_imRK4.F90
                                  source/mod_time_propagator.F90
                                  source/mod_checkpoint.F90
                                  source/mod_imRK4_parameters.F90
                                  source/mod_time_propagator_io.F90
                                  source/mod_cuda.F90
                                  source/mod_nvtx.F90
                                  source/mod_dft.F90
                                  )

if(NOT PLATFORM)
  execute_process(COMMAND bash "-c" "cat /etc/FZJ/systemname 2>/dev/null | tr -d '\n'"
                  OUTPUT_VARIABLE PLATFORM 
                  )
  if(PLATFORM)
    execute_process(COMMAND echo "Setting PLATFORM=${PLATFORM}")
  endif(PLATFORM)
endif()

# Executable names
if(PLATFORM MATCHES jurecadc OR PLATFORM MATCHES juwelsbooster OR PLATFORM MATCHES jusuf)
  set(CMAKE_Fortran_COMPILER mpif90)
elseif(PLATFORM MATCHES osx)
  set(CMAKE_Fortran_COMPILER mpif90)
elseif(PLATFORM MATCHES rwth)
  set(CMAKE_Fortran_COMPILER mpiifort)
elseif(PLATFORM MATCHES jurecabooster)
  set(CMAKE_Fortran_COMPILER mpif90)
elseif(PLATFORM MATCHES juwels)
  set(CMAKE_Fortran_COMPILER mpif90)
elseif(PLATFORM MATCHES iff)
  set(CMAKE_Fortran_COMPILER mpiifort)
elseif(PLATFORM MATCHES hlrs)
  set(CMAKE_Fortran_COMPILER mpif90)
else()
  message( FATAL_ERROR "PLATFORM is not set. Use -DPLATFORM=<machine>" )
endif()

if(NOT PREP)
  set(PREP OFF)
endif()

if(DEBUG)
  set(DEBUGSUFFIX _debug)
  execute_process(COMMAND echo "DEBUG on")
else()
  set(DEBUGSUFFIX)
endif()

if(SUFFIX)
  set(SUFFIX _${SUFFIX})
else()
  set(SUFFIX)
endif()

if(COMPILER MATCHES gfortran OR COMPILER MATCHES gcc OR COMPILER MATCHES gfort OR COMPILER MATCHES gnu)
  set(COMPILER gnu)
  set(CMAKE_Fortran_COMPILER mpif90)
  set(SUFFIX ${SUFFIX}_gnu)
  execute_process(COMMAND echo "Using GCC Compiler")
endif()

if(COMPILER MATCHES psmpi OR COMPILER MATCHES pmpi OR COMPILER MATCHES parastation OR COMPILER MATCHES ps)
  set(SUFFIX ${SUFFIX}_pmpi)
  execute_process(COMMAND echo "Using ParaStationMPI")
endif()

if(COMPILER MATCHES pgi OR COMPILER MATCHES nv OR COMPILER MATCHES nvidia)
  set(COMPILER nv)
  set(SUFFIX ${SUFFIX}_nv)
  execute_process(COMMAND echo "Using NVIDIA Compiler")
  
  set(CMAKE_CXX_COMPILER nvcc)
  set(SOURCE_TITAN_CPP ${SOURCE_TITAN_CPP}  source/diagonalize_gpu.cpp
                                  )
endif()

if(PREP)
  set(CMAKE_Fortran_COMPILER scorep-mpif90)
  set(CMAKE_CXX_COMPILER scorep-nvcc)
  set(INSTRUMENT _inst)
  execute_process(COMMAND echo "Instrumentation on")
else()
  set(INSTRUMENT)
endif()


set(TITAN_EXE titan_${PLATFORM}${DEBUGSUFFIX}${INSTRUMENT}${SUFFIX}.exe)

add_executable(${TITAN_EXE} ${SOURCE_TITAN} ${SOURCE_TITAN_CPP})
set_property(TARGET ${TITAN_EXE} PROPERTY LINKER_LANGUAGE Fortran)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

execute_process(COMMAND git describe --abbrev=8 --dirty --always OUTPUT_VARIABLE GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT COMPILER MATCHES gnu AND NOT COMPILER MATCHES nv)
  if(PLATFORM MATCHES jurecadc OR PLATFORM MATCHES juwelsbooster OR PLATFORM MATCHES jusuf)
    target_link_libraries(${TITAN_EXE} -mkl -qopenmp)
    target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                                -fpp
                                                -DVERSION="${GIT_VERSION}"
                                                -qopenmp
                                                -g
                                                -mkl
                                                -align array64byte
                                                -march=core-avx2
                                                )
  elseif(PLATFORM MATCHES juwels)
    target_link_libraries(${TITAN_EXE} -mkl -qopenmp)
    target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                                -fpp
                                                -xCORE-AVX2
                                                -fma
                                                -qopt-zmm-usage=high
                                                -DVERSION="${GIT_VERSION}"
                                                -qopenmp
                                                -g
                                                -align array64byte
                                                -mkl
                                                )
  elseif(PLATFORM MATCHES osx)
    target_link_libraries(${TITAN_EXE} -mkl -qopenmp)
    target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                                -fpp
                                                -DVERSION="${GIT_VERSION}"
                                                -qopenmp
                                                -g
                                                -mkl
                                                -align array64byte
                                                -xHost
                                                )
  elseif(PLATFORM MATCHES jurecabooster)
    target_link_libraries(${TITAN_EXE} -mkl -qopenmp )
    target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                                -fpp
                                                -DVERSION="${GIT_VERSION}"
                                                -mkl
                                                -g
                                                -xMIC-AVX512
                                                -align array64byte
                                                -qopenmp
                                                )
  elseif(PLATFORM MATCHES rwth)
    target_link_libraries(${TITAN_EXE} -mkl -qopenmp)
    target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                                -fpp
                                                -DVERSION="${GIT_VERSION}"
                                                -qopenmp
                                                -g
                                                -mkl
                                                -align array64byte
                                                -xHost
                                                )
  elseif(PLATFORM MATCHES iff)
    target_link_libraries(${TITAN_EXE} -mkl -qopenmp)
    target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                                -DVERSION="${GIT_VERSION}"
                                                -fpp
                                                -qopenmp
                                                -g
                                                -mkl
                                                -align array64byte
                                                -xHost
                                                )
  elseif(PLATFORM MATCHES hlrs)
    target_link_libraries(${TITAN_EXE} -mkl -qopenmp)
    target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                                -fpp
                                                -DVERSION="${GIT_VERSION}"
                                                -qopenmp
                                                -g
                                                -align array64byte
                                                -mkl
                                                )
  else()
    target_link_libraries(${TITAN_EXE} -mkl -qopenmp)
    target_compile_options(${TITAN_EXE} PRIVATE -no-wrap-margin
                                                -fpp
                                                -DVERSION="${GIT_VERSION}"
                                                -qopenmp
                                                -mkl
                                                #-ipo
                                                #-qopt-report=5
                                                -g
                                                -align array64byte
                                                -xHost
                                                )
  endif(PLATFORM MATCHES jurecadc OR PLATFORM MATCHES juwelsbooster OR PLATFORM MATCHES jusuf)

  if(DEBUG)
    target_compile_definitions(${TITAN_EXE} PRIVATE DEBUG)
    target_compile_options(${TITAN_EXE} PRIVATE -C
                                                -CB
                                                -check uninit
                                                -debug all
                                                -warn all
                                                -e18
                                                -ftrapuv
                                                -traceback
                                                -qopenmp
                                                -fpp
                                                -DVERSION="${GIT_VERSION}"
                                                -g
                                                -O0
                                                )
  else()
    target_compile_options(${TITAN_EXE} PRIVATE -O3)
  endif(DEBUG)

elseif(COMPILER MATCHES gnu) # GFORTRAN:
  target_link_libraries(${TITAN_EXE} -L$ENV{MKLROOT}/lib -Wl,-rpath,$ENV{MKLROOT}/lib -static-libgfortran -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -fopenmp)
  if(NOT DEBUG)
    target_compile_options(${TITAN_EXE} PRIVATE -cpp
                                                -DVERSION="${GIT_VERSION}"
                                                -I$ENV{MKLROOT}/include
                                                -fopenmp
                                                -ffree-form
                                                -ffree-line-length-0
                                                -g
                                                -O3
                                                -march=native
                                                )

  else(NOT DEBUG)
    target_compile_definitions(${TITAN_EXE} PRIVATE DEBUG)
    target_compile_options(${TITAN_EXE} PRIVATE -C
                                                -Wall
                                                -Wextra
                                                -Wconversion
                                                -Wline-truncation
                                                -Wcharacter-truncation
                                                -Wsurprising
                                                -Waliasing
                                                -Wunused-parameter
                                                -fwhole-file
                                                -fimplicit-none
                                                -fbacktrace
                                                -ffree-form
                                                -ffree-line-length-0
                                                -fcheck=all
                                                -ffpe-trap=zero,overflow,underflow
                                                -finit-real=nan
                                                -ftrapv
                                                -fopenmp
                                                -cpp
                                                -DVERSION="${GIT_VERSION}"
                                                -I$ENV{MKLROOT}/include
                                                -g
                                                -O0
                                                )
  endif(NOT DEBUG)

elseif(COMPILER MATCHES nv) # NVIDIA:
  # target_link_libraries(${TITAN_EXE} -L$ENV{MKLROOT}/lib/intel64 -Wl,-rpath,$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -I${CUDA_ROOT}/include -Mcuda -L${CUDA_ROOT}/lib64 -lcusolver -lcudart -lnvToolsExt -acc -gpu=cc70,cc80)

  set_target_properties(${TITAN_EXE} PROPERTIES LINK_FLAGS "-L$ENV{MKLROOT}/lib/intel64 -Wl,-rpath,$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -I${CUDA_ROOT}/include -Mcuda -L${CUDA_ROOT}/lib64 -lcusolver -lcudart -lnvToolsExt -acc -gpu=cc70,cc80")
  if(NOT DEBUG)
    target_compile_definitions(${TITAN_EXE} PRIVATE _GPU)

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
    set_source_files_properties(${SOURCE_TITAN}     PROPERTIES COMPILE_FLAGS ${FORTRAN_FLAGS} )
    string(CONCAT CXX_FLAGS "--gpu-architecture=sm_70,sm_80 "
                            "-lineinfo "
                            "-g "
                            )

    set_source_files_properties(${SOURCE_TITAN_CPP} PROPERTIES COMPILE_FLAGS ${CXX_FLAGS} )

  else(NOT DEBUG)
    target_compile_definitions(${TITAN_EXE} PRIVATE _GPU DEBUG )

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
    set_source_files_properties(${SOURCE_TITAN}     PROPERTIES COMPILE_FLAGS ${FORTRAN_FLAGS})

    string(CONCAT CXX_FLAGS "--gpu-architecture=sm_80 "
                            "-g "
                            "-G "
                            )

    set_source_files_properties(${SOURCE_TITAN_CPP} PROPERTIES COMPILE_FLAGS ${CXX_FLAGS} )

  endif(NOT DEBUG)

endif(NOT COMPILER MATCHES gnu AND NOT COMPILER MATCHES nv)

#target_compile_options(${EXEC_NAME} PRIVATE -no-wrap-margin
#                                        -fpp
#                                        -D _LINUX
#                                        -O3
##                                        -I${NAG}/c_headers#
#         -qopenmp
#         -mkl
#         -xHost
#)

#target_link_libraries(${EXEC_NAME} ${NAG}/lib/libnag_nag.a -qopenmp)

# unset(DEBUG CACHE)
# unset(DEBUGSUFFIX CACHE)
# unset(PREP CACHE)
# unset(SUFFIX CACHE)
# unset(COMPILER CACHE)
