add_subdirectory( reconstruction )

set( vofdir ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof )
set( vof_HEADERS
  brents.hh
  evolution.hh
  flags.hh
  geometricutility.hh
  hyperplane.hh
  reconstruction.hh
  reconstructionSet.hh
  vertexneighborsstencil.hh
)

install(FILES ${vof_HEADERS} DESTINATION ${vofdir})
