cmake_minimum_required(VERSION 3.8)
project(compgeom_hw)

set(CURRENT_TASK A)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../${CURRENT_TASK}/)

set(SOURCE_FILES ${CURRENT_TASK}/main.cpp)
add_executable(${CURRENT_TASK}.o ${SOURCE_FILES})