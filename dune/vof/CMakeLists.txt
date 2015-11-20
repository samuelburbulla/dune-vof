add_subdirectory( reconstruction )

set( vofdir ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof )
set( vof_HEADERS
  brents.hh
  evolution.hh
  femdfwrapper.hh
  flags.hh
  geometricutility.hh
  hyperplane.hh
  mcmgmapper.hh
  reconstruction.hh
  reconstructionSet.hh
  vertexneighborsstencil.hh
)

install(FILES ${vof_HEADERS} DESTINATION ${vofdir})
