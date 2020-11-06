cmake_minimum_required(VERSION 3.12)
project(ToF_example)

enable_language(CXX)

#set(CMAKE_CXX_STANDARD 14)

find_package(MythiCaL REQUIRED)
add_executable(ToF_example ToF_example.cpp)

target_link_libraries(ToF_example PUBLIC mythical::mythical)
