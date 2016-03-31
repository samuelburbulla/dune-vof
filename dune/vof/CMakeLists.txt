add_subdirectory( geometry )
add_subdirectory( reconstruction )
add_subdirectory( stencil )

set( vofdir ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof )
set( vof_HEADERS
  brents.hh
  evolution.hh
  femdfwrapper.hh
  flags.hh
  mcmgmapper.hh
  reconstruction.hh
  reconstructionSet.hh
)

install(FILES ${vof_HEADERS} DESTINATION ${vofdir})
