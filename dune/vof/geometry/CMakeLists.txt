set( geometrydir ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof/geometry )
set( geometry_HEADERS
  2d/polygon.hh
  2d/intersect.hh
  algorithm.hh
  halfspace.hh
  intersect.hh
  utility.hh
)

install(FILES ${geometry_HEADERS} DESTINATION ${geometrydir})
