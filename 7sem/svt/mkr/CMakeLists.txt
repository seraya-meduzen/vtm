cmake_minimum_required(VERSION 3.20)
project(svt)

find_package(inmost REQUIRED)

link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})

add_executable(svt main.cpp)
target_link_libraries(svt ${INMOST_LIBRARIES})