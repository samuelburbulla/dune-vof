set(HEADERS
  heightfunction.hh
  generalheightfunction.hh
  swartz.hh
  modifiedswartz.hh
  modifiedyoungs.hh
  modifiedyoungssecondorder.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof/reconstruction)
