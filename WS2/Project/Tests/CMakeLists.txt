cmake_minimum_required(VERSION 3.9)
project(tests)

find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(CMAKE_C_STANDARD 99)

set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wall -pedantic -Werror -Wno-write-strings")

# Set the normal target
set(TEST_SOURCE_FILES testing.c uvComm.test.c getProcessCoordinates.test.c getProcessNeighbours.test.c pressureComm.test.c)
add_executable(tests tests.c ${TEST_SOURCE_FILES})
target_link_libraries(tests core ${MPI_C_LIBRARIES} m)

# The below is to always get an updated copy of cavity100.dat inside the cmake-build-debug folder where the binary is.
#add_custom_target(copy_aux_files COMMAND cp *.dat *.pgm ${sim_BINARY_DIR}/ WORKING_DIRECTORY ${sim_SOURCE_DIR})
#add_dependencies(sim copy_aux_files)
#add_custom_command(TARGET sim POST_BUILD COMMAND cp cavity100.dat ${sim_BINARY_DIR}/ WORKING_DIRECTORY ${sim_SOURCE_DIR})
#add_custom_command(OUTPUT execute_always COMMAND cp cavity100.dat ${sim_BINARY_DIR}/
#        WORKING_DIRECTORY ${sim_SOURCE_DIR} DEPENDS ${sim_SOURCE_DIR}/cavity100.dat)
