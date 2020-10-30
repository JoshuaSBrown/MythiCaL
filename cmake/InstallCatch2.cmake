# If Catch2 is not found, instead use the git submodule
if (NOT EXISTS ${mythical_SOURCE_DIR}/external/Catch2/single_include)
  # Unable to find the header files for Catch2 or they don't exist
  message(STATUS "Downloading Catch2 submodule.")

  # Clone the submodule
  execute_process(COMMAND git submodule update --init --force -- external/Catch2 WORKING_DIRECTORY ${mythical_SOURCE_DIR})
endif()

add_subdirectory(${mythical_SOURCE_DIR}/external/Catch2)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/external/Catch2/contrib")
