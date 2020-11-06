cmake_minimum_required(VERSION 3.12)
project(print_version)

find_package(MythiCaL REQUIRED)
add_executable(example example.cpp)

target_link_libraries(example PUBLIC mythical::mythical)
