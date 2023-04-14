cmake_minimum_required(VERSION 3.22)
PROJECT (Hello)

# set(CMAKE_CXX_COMPILER mpic++)
set(CMAKE_CXX_COMPILER "/usr/bin/mpic++")
set(GCC_COVERAGE_COMPILE_FLAGS "-flto -Ofast -march=native -llapack -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS} -lblas")


find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})


add_executable(hello toom-cook-MPI.cpp)
target_link_libraries(hello ${MPI_LIBRARIES})

# target_link_libraries(hello OpenMP::OpenMP_CXX)
