set( reconstructiondir ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof/reconstruction )
set( reconstruction_HEADERS
  modifiedswartz.hh
  modifiedyoungs.hh
)

install(FILES ${vof_HEADERS} DESTINATION ${reconstructiondir})
