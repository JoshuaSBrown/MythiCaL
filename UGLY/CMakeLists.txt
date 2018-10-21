cmake_minimum_required (VERSION 3.2)

enable_language(CXX)

project (ugly)

########################################################
# Compiler Flags                                       #
########################################################

IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -pedantic")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

########################################################
# Basic system tests (standard libraries, headers etc) #
########################################################
set(CMAKE_REQUIRED_INCLUDES "/usr/local/include")

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

###############################
# Add header files to project #
###############################

file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/ugly/*.hpp)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/libugly/edge)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/libugly/graphvisitor)

# Only install the public facing header files
install(FILES ${HEADERS} DESTINATION include/ugly)

###############################
# Add source files to project #
###############################

file(GLOB SOURCES1 ${CMAKE_CURRENT_SOURCE_DIR}/src/libugly/*.cpp)
file(GLOB SOURCES2 ${CMAKE_CURRENT_SOURCE_DIR}/src/libugly/edge/*.cpp)
file(GLOB SOURCES3 ${CMAKE_CURRENT_SOURCE_DIR}/src/libugly/edge/*.hpp)
file(GLOB SOURCES4 ${CMAKE_CURRENT_SOURCE_DIR}/src/libugly/graphvisitor/*.cpp)
file(GLOB SOURCES5 ${CMAKE_CURRENT_SOURCE_DIR}/src/libugly/graphvisitor/*.hpp)
add_library(ugly ${SOURCES1} ${SOURCES2} ${SOURCES3} ${SOURCES4} ${SOURCES5})
install(TARGETS ${PROJECT_NAME} DESTINATION lib/${PROJECT_NAME})

set_target_properties(ugly PROPERTIES LINKER_LANGUAGE CXX)

##########################
# Check if testing is on #
##########################
option(BUILD_SHARED_LIBS "Build shared libs" ON)
option(ENABLE_TESTING "Build and enable testing stuff" OFF)
if(ENABLE_TESTING)
  enable_testing()
  add_subdirectory(src/tests)
endif(ENABLE_TESTING)
