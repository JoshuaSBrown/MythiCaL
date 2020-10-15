cmake_minimum_required(VERSION 3.4)

enable_language(CXX)

project(mythical)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_SHARED_LIBS "Build shared libs" ON)
option(ENABLE_TESTING "Build and enable unit testing" OFF)
########################################################
# Compiler Flags                                       #
########################################################

IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -pedantic")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)
IF(CMAKE_BUILD_TYPE MATCHES "DEBUG")
  message("Building in Debug mode")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
ENDIF(CMAKE_BUILD_TYPE MATCHES "DEBUG")

###########################
# Include Extra Functions #
###########################
include(GNUInstallDirs)
include(CheckIncludeFileCXX)
include(CheckLibraryExists)
########################################################
# Basic system tests (standard libraries, headers etc) #
########################################################

foreach(HEADER cmath cassert fstream functional iostream limits list map ostream sstream stack stdexcept string vector)
  check_include_file_cxx(${HEADER} FOUND_${HEADER})
  if(NOT FOUND_${HEADER})
    message(FATAL_ERROR "Could not find needed header - ${HEADER}")
  endif(NOT FOUND_${HEADER})
endforeach(HEADER)

set(MATH_LIBRARIES "m" CACHE STRING "math library")
mark_as_advanced( MATH_LIBRARIES )
foreach(FUNC sqrt)
  check_library_exists(${MATH_LIBRARIES} ${FUNC} "" FOUND_${FUNC}_${MATH_LIBRARIES})
  if(NOT FOUND_${FUNC}_${MATH_LIBRARIES})
    message(FATAL_ERROR "Could not find needed math function - ${FUNC}")
  endif(NOT FOUND_${FUNC}_${MATH_LIBRARIES})
endforeach(FUNC)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/CoarseGrainSites/src/libmythical/topologyfeatures ${CMAKE_CURRENT_SOURCE_DIR}/CoarseGrainSites/include )

###########################
# Setup relative includes #
###########################
include_directories(CoarseGrainSites/include)
include_directories(CoarseGrainSites/src)
include_directories(UGLY/include)

########################################################################
# Grab source files and set the directory where they will be installed #
########################################################################
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/CoarseGrainSites/src/libmythical/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/CoarseGrainSites/src/libmythical/topologyfeatures/*.cpp)
file(GLOB SOURCES_UGLY1 ${CMAKE_CURRENT_SOURCE_DIR}/UGLY/src/libugly/*.cpp)
file(GLOB SOURCES_UGLY2 ${CMAKE_CURRENT_SOURCE_DIR}/UGLY/src/libugly/edge/*.cpp)
file(GLOB SOURCES_UGLY3 ${CMAKE_CURRENT_SOURCE_DIR}/UGLY/src/libugly/edge/*.hpp)
file(GLOB SOURCES_UGLY4 ${CMAKE_CURRENT_SOURCE_DIR}/UGLY/src/libugly/graphvisitor/*.cpp)
file(GLOB SOURCES_UGLY5 ${CMAKE_CURRENT_SOURCE_DIR}/UGLY/src/libugly/graphvisitor/*.hpp)
add_library(mythical ${SOURCES} ${SOURCES_UGLY1} ${SOURCES_UGLY2} ${SOURCES_UGLY3} ${SOURCES_UGLY4} ${SOURCES_UGLY5})
set_target_properties(mythical PROPERTIES LINKER_LANGUAGE CXX)

add_subdirectory(CoarseGrainSites/include/mythical)
include(cmake/MythiCaLInstall.cmake)
###############################
# Check if unit testing is on #
###############################
if(ENABLE_TESTING)
  enable_testing()
  message("testing enabled")
  IF(CMAKE_BUILD_TYPE MATCHES "RELEASE")
    message("Cannot build unit tests with release mode")  
  ELSE(CMAKE_BUILD_TYPE MATCHES "RELEASE")
    add_subdirectory(CoarseGrainSites/src/tests/unit)
  ENDIF(CMAKE_BUILD_TYPE MATCHES "RELEASE")
  add_subdirectory(CoarseGrainSites/src/tests/regression)
endif(ENABLE_TESTING)

###########################################
# Check if performance testing is enabled #
###########################################

if(CXXTEST_ADD_PERFORMANCE)
  enable_testing()
  message("Performance testing enabled")
  add_subdirectory(CoarseGrainSites/src/tests/performance)
endif(CXXTEST_ADD_PERFORMANCE)

#####################################
# Check if code coverage is enabled #
#####################################

foreach(FLAG_OPT ${CMAKE_CXX_FLAGS})
  if(${FLAG_OPT} MATCHES "-coverage")
    find_program(GCOV_EXE gcov)
    if(GCOV_EXE)
      file(GLOB COVERAGE_FILES ${PROJECT_BINARY_DIR}/CMakeFiles/mythical.dir/CoarseGrainSites/src/libmythical/*gcno ${PROJECT_BINARY_DIR}/CMakeFiles/mythical.dir/CoarseGrainSites/src/libmythical/topologyfeatures/*gcno)
      foreach(PROG ${COVERAGE_FILES})
        get_filename_component(NAME_FILE ${PROG} NAME)
        add_test(coverage_${NAME_FILE} gcov ${PROG})
      endforeach()

    elseif(NOT GCOV_EXE)
      message("Unable to find gcov")
    endif(GCOV_EXE)
  endif(${FLAG_OPT} MATCHES "-coverage")
endforeach()
