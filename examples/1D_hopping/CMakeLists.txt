cmake_minimum_required(VERSION 3.12)
project(1D_example)

enable_language(CXX)

find_package(MythiCaL REQUIRED)
add_executable(1D_example 1D_example.cpp)

target_link_libraries(1D_example PUBLIC mythical::mythical)
