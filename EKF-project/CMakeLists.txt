cmake_minimum_required(VERSION 3.10.2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(get_telemetry)

add_executable(get_telemetry
    get_telemetry.cpp
)

find_package(MAVSDK REQUIRED)

target_link_libraries(get_telemetry
    MAVSDK::mavsdk
)

if(NOT MSVC)
    add_compile_options(get_telemetry PRIVATE -Wall -Wextra)
else()
    add_compile_options(get_telemetry PRIVATE -WX -W2)
endif()
