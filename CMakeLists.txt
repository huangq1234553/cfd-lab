cmake_minimum_required(VERSION 3.9)
project(sim)

set(CMAKE_C_STANDARD 99)
#set(PRECICE_ROOT /home/keefe/precice-1.1.1)
#set(PRECICE_ROOT ENV{PRECICE_ROOT})
set(PRECICE_ROOT /opt/preCICE/precice-1.1.1)
set(INCLUDES -I${PRECICE_ROOT}/src/precice)
set(LIBS precice)
set(GCC_LINKER_FLAGS "-L${PRECICE_ROOT}/build/last")

set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wall -pedantic -Werror ${INCLUDES}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_LINKER_FLAGS}")

set(SOURCE_FILES main.c boundary_val.c helper.c init.c sor.c uvp.c visual.c logger.c boundary_configurator.c timing.c precice_adapter.c)
add_executable(sim ${SOURCE_FILES})
target_link_libraries(sim ${LIBS} m)
include_directories(${PRECICE_ROOT}/src/precice)


# The below is to always get an updated copy of cavity100.dat inside the cmake-build-debug folder where the binary is.
#add_custom_target(copy_aux_files COMMAND cp *.dat *.pgm ${sim_BINARY_DIR}/ WORKING_DIRECTORY ${sim_SOURCE_DIR})
#add_dependencies(sim copy_aux_files)
##add_custom_command(TARGET sim POST_BUILD COMMAND cp cavity100.dat ${sim_BINARY_DIR}/ WORKING_DIRECTORY ${sim_SOURCE_DIR})
##add_custom_command(OUTPUT execute_always COMMAND cp cavity100.dat ${sim_BINARY_DIR}/
#        WORKING_DIRECTORY ${sim_SOURCE_DIR} DEPENDS ${sim_SOURCE_DIR}/cavity100.dat)
