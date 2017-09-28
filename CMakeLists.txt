cmake_minimum_required(VERSION 3.8)
project(compgeom_hw)

set(CURRENT_TASK c)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../${CURRENT_TASK}/)

file(GLOB SOURCE_FILES ${CURRENT_TASK}/*.cpp)
include_directories(headers/)

add_executable(main ${SOURCE_FILES})