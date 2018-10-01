cmake_minimum_required(VERSION 3.3)

enable_language(CXX)

project(kmc_course_grain)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

########################################################
# Basic system tests (standard libraries, headers etc) #
########################################################

include(CheckIncludeFile)
foreach(HEADER assert.h math.h stdio.h stdlib.h string.h)
  check_include_file(${HEADER} FOUND_${HEADER})
  if(NOT FOUND_${HEADER})
    message(FATAL_ERROR "Could not find needed header - ${HEADER}")
  endif(NOT FOUND_${HEADER})
endforeach(HEADER)

include(CheckIncludeFileCXX)
foreach(HEADER cmath fstream functional iostream limits list map ostream sstream stack stdexcept string vector)
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

file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/kmcoursegrain/*.hpp)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
install(FILES ${HEADERS} DESTINATION include/kmccoursegrain)

file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/libkmccoursegrain/*.cpp)
add_library(kmc_course_grain ${SOURCES})
install(TARGETS ${PROJECT_NAME} DESTINATION lib/${PROJECT_NAME})

set_target_properties(kmc_course_grain PROPERTIES LINKER_LANGUAGE CXX)

##########################
# Check if testing is on #
##########################
option(BUILD_SHARED_LIBS "Build shared libs" ON)
option(ENABLE_TESTING "Build and enable testing stuff" OFF)
if(ENABLE_TESTING)
  enable_testing()
  add_subdirectory(src/tests)
endif(ENABLE_TESTING)

