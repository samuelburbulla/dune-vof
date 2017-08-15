add_subdirectory(common)
add_subdirectory(curvature)
add_subdirectory(evolution)
add_subdirectory(geometry)
add_subdirectory(interfacegrid)
add_subdirectory(interpolation)
add_subdirectory(reconstruction)
add_subdirectory(stencil)
add_subdirectory(test)

set(HEADERS
  algorithm.hh
  brents.hh
  colorfunction.hh
  curvatureset.hh
  dataset.hh
  eoc.hh
  evolution.hh
  flagging.hh
  flagset.hh
  mixedcellmapper.hh
  reconstruction.hh
  reconstructionset.hh
  utility.hh
  velocity.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof)
