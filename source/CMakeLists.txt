
file(GLOB sourceFiles RELATIVE "${CMAKE_PROJECT_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/*.F90")
#message(STATUS "${sourceFiles}")

set(SOURCE_TITAN ${SOURCE_TITAN} ${sourceFiles} PARENT_SCOPE )
