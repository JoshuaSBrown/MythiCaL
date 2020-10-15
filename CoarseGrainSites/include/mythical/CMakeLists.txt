# Maintain directory structure in installed include files
install(DIRECTORY . DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/mythical" FILES_MATCHING PATTERN "*.hpp")
