cmake_minimum_required(VERSION 3.10)

# Project Information
project(PROJECT_HEAT_CYLINDER3D)

# Set MPI as the compiler
find_package(MPI REQUIRED)
if (MPI_FOUND)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

# C/C++ Standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include Directory
include_directories(${PROJECT_SOURCE_DIR}/include)

# Add O3 optimization flag
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Add OpenMP flag
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# Executable Files
add_executable(main main.cc)
add_executable(cylinder main_cylinder.cc)