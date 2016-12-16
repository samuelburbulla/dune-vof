set(HEADERS
  edgeneighborsstencil.hh
  vertexneighborsstencil.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof/stencil)
