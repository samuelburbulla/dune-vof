configure_file( parameter.ini ${CMAKE_CURRENT_BINARY_DIR}/parameter.ini COPYONLY )
configure_file( test8.dgf ${CMAKE_CURRENT_BINARY_DIR}/test8.dgf COPYONLY )

add_executable("dune-vof" main.cc polygon.hh vtu.hh)
target_link_dune_default_libraries("dune-vof")
set_property( TARGET dune-vof APPEND PROPERTY
  COMPILE_DEFINITIONS "ALUGRID_CUBE;GRIDDIM=2" )
