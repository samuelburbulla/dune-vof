set( vofdir ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof )
set( vof_HEADERS
  ellipse.hh
  linearwall.hh
  problems.hh
  rotatingcircle.hh
  sflow.hh
  slottedcylinder.hh
)

install(FILES ${vof_HEADERS} DESTINATION ${vofdir})
