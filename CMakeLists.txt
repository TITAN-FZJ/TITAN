# Projects name
project(TITAN)

# Minimum required Version
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

message( STATUS ">>-----------------------------------------------<<")
message( STATUS ">>---------------------TITAN---------------------<<")

option(DEBUG "Enable debugging" OFF)

enable_language(Fortran)

set(CMAKE_Fortran_COMPILER mpiifort)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(NAG /opt/NAG/fll6i26dcl)
set(SOURCE_TITAN)
set(EXEC_NAME main.exe)
if(DEBUG)
  set(EXEC_NAME main_debug.exe)
endif(DEBUG)

include_directories(${PROJECT_SOURCE_DIR}/source)

add_subdirectory(${PROJECT_SOURCE_DIR}/source)


add_executable(main.exe ${SOURCE_TITAN})

target_compile_options(main.exe PRIVATE -no-wrap-margin
                                        -fpp
                                        -D _LINUX
                                        -g
                                        -O0
                                        -warn all
                                        -I${NAG}/c_headers)
if(DEBUG)
target_compile_options(${EXEC_NAME} PRIVATE -CB
                                            -check-all
                                            -check uninit
                                            -ftrapuv
                                            -debug all
                                            -traceback
                                            -g
                                            -warn all
                                            -O0)
endif(DEBUG)

target_link_libraries(main.exe ${NAG}/lib/libnag_nag.a)
