cmake_minimum_required(VERSION 3.4)

enable_language(CXX)

project(kmc_course_grain)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
########################################################
# Compiler Flags                                       #
########################################################

IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -pedantic")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

########################################################
# Add and build other dependent libraries first        #
########################################################

add_subdirectory(UGLY)
########################################################
# Basic system tests (standard libraries, headers etc) #
########################################################

include(CheckIncludeFileCXX)
foreach(HEADER cmath cassert fstream functional iostream limits list map ostream sstream stack stdexcept string vector)
  check_include_file_cxx(${HEADER} FOUND_${HEADER})
  if(NOT FOUND_${HEADER})
    message(FATAL_ERROR "Could not find needed header - ${HEADER}")
  endif(NOT FOUND_${HEADER})
endforeach(HEADER)

set(MATH_LIBRARIES "m" CACHE STRING "math library")
mark_as_advanced( MATH_LIBRARIES )
include(CheckLibraryExists)
foreach(FUNC sqrt)
  check_library_exists(${MATH_LIBRARIES} ${FUNC} "" FOUND_${FUNC}_${MATH_LIBRARIES})
  if(NOT FOUND_${FUNC}_${MATH_LIBRARIES})
    message(FATAL_ERROR "Could not find needed math function - ${FUNC}")
  endif(NOT FOUND_${FUNC}_${MATH_LIBRARIES})
endforeach(FUNC)

file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/CourseGrainSites/include/kmcoursegrain/*.hpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/CourseGrainSites/src/libkmccoursegrain/topologyfeatures ${CMAKE_CURRENT_SOURCE_DIR}/CourseGrainSites/include )

install(FILES ${HEADERS} DESTINATION CourseGrainSites/include/kmccoursegrain)

message(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/CourseGrainSites/src/libkmccoursegrain/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/CourseGrainSites/src/libkmccoursegrain/topologyfeatures/*.cpp)

add_library(kmc_course_grain ${SOURCES} )
install(TARGETS ${PROJECT_NAME} DESTINATION lib/${PROJECT_NAME})

set_target_properties(kmc_course_grain PROPERTIES LINKER_LANGUAGE CXX)

##########################
# Check if testing is on #
##########################
option(BUILD_SHARED_LIBS "Build shared libs" ON)

option(ENABLE_TESTING "Build and enable unit testing" OFF)
if(ENABLE_TESTING)
  enable_testing()
  message("Unit testing enabled")
  add_subdirectory(CourseGrainSites/src/tests/unit)
endif(ENABLE_TESTING)

if(CXXTEST_ADD_PERFORMANCE)
  enable_testing()
  message("Performance testing enabled")
  add_subdirectory(CourseGrainSites/src/tests/performance)
endif(CXXTEST_ADD_PERFORMANCE)
