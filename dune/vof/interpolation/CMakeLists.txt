set(HEADERS
  interpolations.hh
  circleinterpolation.hh
  recursiveinterpolation.hh
  recursiveinterpolationcube.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof/interpolation)
