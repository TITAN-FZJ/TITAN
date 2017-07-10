# Minimum required Version
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

# Projects name
project(TITAN)

option(DEBUG "Enable debugging" OFF)

enable_language(Fortran)

file(GLOB sources source/*.F90
	  	  source/*.f90)

add_executable(titan.exe ${sources})

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
                                            -check-all
                                            -debug all
                                            -warn all
                                            -ftrapuv
                                            -traceback
                                            -g
                                            -O0
                      )
endif(DEBUG)

#target_link_libraries(${EXEC_NAME} ${NAG}/lib/libnag_nag.a -qopenmp)
