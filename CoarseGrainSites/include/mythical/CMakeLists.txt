configure_file(version.hpp.in version.hpp @ONLY)
# Maintain directory structure in installed include files
install(DIRECTORY . DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/mythical" FILES_MATCHING PATTERN "*.hpp")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.hpp DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/mythical")
