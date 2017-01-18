add_subdirectory(curvature)
add_subdirectory(geometry)
add_subdirectory(reconstruction)
add_subdirectory(stencil)
add_subdirectory(test)

set(HEADERS
  brents.hh
  curvatureSet.hh
  dataset.hh
  evolution.hh
  interpolation.hh
  mixedcellmapper.hh
  flags.hh
  reconstruction.hh
  reconstructionSet.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof)
