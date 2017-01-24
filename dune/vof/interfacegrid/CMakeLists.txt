set(HEADERS
  capabilities.hh
  datahandle.hh
  declaration.hh
  entity.hh
  entityseed.hh
  geometry.hh
  grid.hh
  gridview.hh
  hostgridaccess.hh
  idset.hh
  indexset.hh
  intersection.hh
  intersectioniterator.hh
  iterator.hh
  persistentcontainer.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof/interfacegrid)
