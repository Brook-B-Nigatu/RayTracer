cmake_minimum_required(VERSION 3.10)

project(raytracer)

set(CMAKE_CXX_STANDARD 11)

add_executable(main TD1/ray_tracer.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(main PUBLIC OpenMP::OpenMP_CXX)
endif()

target_include_directories(main PRIVATE include)