cmake_minimum_required(VERSION 3.22)
project(svtfortran)
enable_language(Fortran)

set(GCC_COVERAGE_COMPILE_FLAGS "-Ofast -march=native")
set(GCC_COVERAGE_COMPILE_FLAGS "-g -fsanitize=address")
set(CMAKE_Fortran_FLAGS "${GCC_COVERAGE_COMPILE_FLAGS} -lblas")

find_package(BLAS REQUIRED)

add_executable(svtspars main.f90)
target_link_libraries(svtspars ${CMAKE_SOURCE_DIR}/libskit.a)
