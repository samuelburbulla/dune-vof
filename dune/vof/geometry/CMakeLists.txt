add_subdirectory(2d)
add_subdirectory(3d)

set(HEADERS
  algorithm.hh
  halfspace.hh
  intersect.hh
  polytope.hh
  upwindpolygon.hh
  utility.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof/geometry)
