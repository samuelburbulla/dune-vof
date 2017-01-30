add_subdirectory(curvature)
add_subdirectory(geometry)
add_subdirectory(interfacegrid)
add_subdirectory(reconstruction)
add_subdirectory(stencil)
add_subdirectory(test)

set(HEADERS
  brents.hh
  curvatureset.hh
  dataset.hh
  evolution.hh
  interpolation.hh
  mixedcellmapper.hh
  reconstruction.hh
  reconstructionset.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof)
