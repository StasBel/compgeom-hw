cmake_minimum_required(VERSION 3.8)
project(compgeom_hw)

set(CURRENT_TASK a)

set(CMAKE_C_COMPILER /usr/local/Cellar/gcc@5/5.4.0_1/bin/gcc-5 CACHE PATH "" FORCE)
set(CMAKE_CXX_COMPILER /usr/local/Cellar/gcc@5/5.4.0_1/bin/g++-5 CACHE PATH "" FORCE)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../${CURRENT_TASK}/)

file(GLOB SOURCE_FILES ${CURRENT_TASK}/*.cpp)
add_executable(main ${SOURCE_FILES})