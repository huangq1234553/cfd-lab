cmake_minimum_required(VERSION 3.9)
project(sim)

set(CMAKE_C_STANDARD 99)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Werror")

set(SOURCE_FILES main.c boundary_val.c helper.c init.c sor.c uvp.c visual.c)
add_executable(sim ${SOURCE_FILES})
target_link_libraries(sim m)
add_custom_command(TARGET sim POST_BUILD COMMAND cp cavity100.dat ${sim_BINARY_DIR}/ WORKING_DIRECTORY ${sim_SOURCE_DIR})
#add_custom_target(sim COMMAND make -C ${sim_SOURCE_DIR})
