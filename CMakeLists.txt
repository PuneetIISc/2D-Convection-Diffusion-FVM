cmake_minimum_required(VERSION 3.10)
project(cd2d_unstructured CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(cd2d
    src/main.cpp
    src/Input.cpp
    src/Mesh.cpp
    src/Solver.cpp
    src/TimeIntegrator.cpp
    src/Writer.cpp
)
target_include_directories(cd2d PRIVATE include)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(cd2d PRIVATE -O2 -Wall)
endif()

configure_file(${CMAKE_SOURCE_DIR}/userchoice.inp ${CMAKE_BINARY_DIR}/userchoice.inp COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/file.input      ${CMAKE_BINARY_DIR}/file.input      COPYONLY)
