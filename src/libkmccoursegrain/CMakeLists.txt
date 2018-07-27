
#file(GLOB KMCCOURSEGRAIN_SOURCES *.cpp)
#
#include_directories(${CMAKE_CURRENT_BINARY_DIR})
#add_library(kmc_course_grain SHARED ${KMCCOURSEGRAIN_SOURCES})
#
#INSTALL(
#  TARGETS kmc_course_grain
#  ARCHIVE DESTINATION lib
#  LIBRARY DESTINATION lib  # <-- Add this line
#  COMPONENT library
#)
#
#set_target_properties(kmc_course_grain PROPERTIES SOVERSION ${SOVERSION})
#install(TARGETS kmc_course_grain LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
#
