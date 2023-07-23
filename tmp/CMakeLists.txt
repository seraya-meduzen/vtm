cmake_minimum_required(VERSION 3.22)
PROJECT (Hello)

set(GCC_COVERAGE_COMPILE_FLAGS "-flto -Ofast -march=native -llapack -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS} -lblas")


find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# message("${LAPACK_LIBRIRIES}")
# message("${LAPACK_FOUND}")


add_executable(hello toom-cook.cpp)
target_link_libraries(hello OpenMP::OpenMP_CXX)
