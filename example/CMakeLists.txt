set( vof_HEADERS
  average.hh
  colorFunction.hh
  errors.hh
  io.hh
  polygon.hh
  problem.hh
  vtu.hh
)

configure_file( parameter.ini ${CMAKE_CURRENT_BINARY_DIR}/parameter.ini COPYONLY )
configure_file( cube2d.dgf ${CMAKE_CURRENT_BINARY_DIR}/cube2d.dgf COPYONLY )

add_executable( vof main.cc)
target_link_dune_default_libraries(vof)
set_property( TARGET vof APPEND PROPERTY
  COMPILE_DEFINITIONS "ALUGRID_SIMPLEX;GRIDDIM=2" )
