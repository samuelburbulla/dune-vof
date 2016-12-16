set(HEADERS
  intersect.hh
  polygon.hh
  upwindpolygon.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof/geometry/2d)
