set( stencildir ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof/stencil )
set( stencil_HEADERS
  edgeneighborsstencil.hh
  vertexneighborsstencil.hh
)

install(FILES ${stencil_HEADERS} DESTINATION ${stencildir})
