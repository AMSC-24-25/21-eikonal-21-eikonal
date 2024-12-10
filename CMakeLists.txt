cmake_minimum_required(VERSION 3.10)
project(EikonalSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_path(EIGEN3_INCLUDE_DIR Eigen/Eigen
  HINTS ${mkEigenInc} $ENV{mkEigenInc})

include_directories(${EIGEN3_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include)

add_definitions(-DDIMENSION=2)

file(GLOB LOCAL_PROBLEM_SOURCES "LocalProblem/*.cpp")
list(REMOVE_ITEM LOCAL_PROBLEM_SOURCES "${CMAKE_SOURCE_DIR}/LocalProblem/main_eikonal.cpp")

add_executable(main_global_solver src/main_global_solver.cpp ${LOCAL_PROBLEM_SOURCES})
