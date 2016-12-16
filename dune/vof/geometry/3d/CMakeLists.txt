set(HEADERS
  face.hh
  intersect.hh
  polygonwithdirections.hh
  polyhedron.hh
  rotation.hh
  upwindpolygon.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof/geometry/3d)
