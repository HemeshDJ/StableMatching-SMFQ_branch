cmake_minimum_required(VERSION 3.0)
project(graphmatching)

# compiler and flags
#add_definitions("-Wall -O3")
add_definitions("-Wall -O0 -g")
# -fno-limit-debug-info")

set(CMAKE_CXX_STANDARD 17)

# include directory and the libraries we will build
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
add_subdirectory(lib)
# add_subdirectory(test)

# build
add_executable(graphmatching main.cc)
target_link_libraries(graphmatching graph)

add_executable(tester test.cc)
target_link_libraries(tester graph)

install(TARGETS graphmatching RUNTIME DESTINATION bin)
install(TARGETS tester RUNTIME DESTINATION bin)