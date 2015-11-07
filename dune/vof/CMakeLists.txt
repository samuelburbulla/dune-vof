set( vofdir ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof )
set( vof_HEADERS
  brents.hh
  evolution.hh
  flags.hh
  geometricutility.hh
  hyperplane.hh
  modifiedswartz.hh
  reconstruct.hh
  reconstructionSet.hh
  secondyoungsnormalguessing.hh
  vertexneighborsstencil.hh
)

install(FILES ${vof_HEADERS} DESTINATION ${vofdir})
