add_subdirectory(curvature)
add_subdirectory(geometry)
add_subdirectory(interfacegrid)
add_subdirectory(reconstruction)
add_subdirectory(stencil)
add_subdirectory(test)

set(HEADERS
  algorithm.hh
  brents.hh
  curvatureset.hh
  dataset.hh
  evolution.hh
  flagging.hh
  flagset.hh
  interpolation.hh
  mixedcellmapper.hh
  reconstruction.hh
  reconstructionset.hh
  utility.hh
  velocity.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof)
