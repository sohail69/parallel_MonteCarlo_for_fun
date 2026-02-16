# CMake minimum version
cmake_minimum_required (VERSION 3.20)

# Project name
project(pFEM
        DESCRIPTION "A Hierachical structured grid Monte-Carlo sampling framework"
        LANGUAGES C CXX)

add_executable(monteCarlo mainMonteCarloTest.cpp)

# Import MPI and OpenMP
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(${MPI_INCLUDE_PATH})

# Generate obj files
target_link_libraries(monteCarlo PRIVATE OpenMP::OpenMP_CXX ${MPI_LIBRARIES})

