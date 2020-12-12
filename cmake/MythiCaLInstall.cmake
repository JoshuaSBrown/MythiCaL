
set_target_properties(mythical PROPERTIES VERSION ${mythical_VERSION} SOVERSION ${SOVERSION})

target_include_directories(mythical PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/libmythical>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/mythical>
  )

message("Instal include dir ${CMAKE_INSTALL_INCLUDEDIR}/mythical")

install(TARGETS mythical EXPORT mythicalTargets
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  )

# Install config files
install(FILES ${PROJECT_SOURCE_DIR}/cmake/MythiCaLConfig.cmake DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/mythical")

install(EXPORT mythicalTargets
  FILE mythicalTargets.cmake
  NAMESPACE mythical::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/mythical"
  )
