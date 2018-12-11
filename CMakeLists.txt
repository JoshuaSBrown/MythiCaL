cmake_minimum_required(VERSION 3.4)

enable_language(CXX)

project(kmc_coarse_grain)

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

file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/CoarseGrainSites/include/kmccoarsegrain/*.hpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/CoarseGrainSites/src/libkmccoarsegrain/topologyfeatures ${CMAKE_CURRENT_SOURCE_DIR}/CoarseGrainSites/include )

install(FILES ${HEADERS} DESTINATION include/kmccoarsegrain)

message( ${CMAKE_CURRENT_SOURCE_DIR} )

file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/CoarseGrainSites/src/libkmccoarsegrain/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/CoarseGrainSites/src/libkmccoarsegrain/topologyfeatures/*.cpp)

add_library(kmc_coarse_grain ${SOURCES} )
install(TARGETS ${PROJECT_NAME} DESTINATION lib/${PROJECT_NAME})

set_target_properties(kmc_coarse_grain PROPERTIES LINKER_LANGUAGE CXX)

###############################
# Check if unit testing is on #
###############################
option(BUILD_SHARED_LIBS "Build shared libs" ON)

option(ENABLE_TESTING "Build and enable unit testing" OFF)
if(ENABLE_TESTING)
  enable_testing()
  message("Unit testing enabled")
  add_subdirectory(CoarseGrainSites/src/tests/unit)
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
      file(GLOB COVERAGE_FILES ${PROJECT_BINARY_DIR}/CMakeFiles/kmc_coarse_grain.dir/CoarseGrainSites/src/libkmccoarsegrain/*gcno ${PROJECT_BINARY_DIR}/CMakeFiles/kmc_coarse_grain.dir/CoarseGrainSites/src/libkmccoarsegrain/topologyfeatures/*gcno)
      foreach(PROG ${COVERAGE_FILES})
        get_filename_component(NAME_FILE ${PROG} NAME)
        add_test(coverage_${NAME_FILE} gcov ${PROG})
      endforeach()

    elseif(NOT GCOV_EXE)
      message("Unable to find gcov")
    endif(GCOV_EXE)
  endif(${FLAG_OPT} MATCHES "-coverage")
endforeach()
