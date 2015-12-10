
# copy data to build source to make tests work
configure_file( 1dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/1dgrid.dgf COPYONLY )
configure_file( 2dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY )
configure_file( 3dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/3dgrid.dgf COPYONLY )
configure_file( parameter ${CMAKE_CURRENT_BINARY_DIR}/parameter COPYONLY )


add_executable( "vof-fem" main.cc ${HEADERS} )
target_link_dune_default_libraries( "vof-fem" )

set_property( TARGET vof-fem APPEND PROPERTY
  COMPILE_DEFINITIONS "ALUGRID_CUBE;GRIDDIM=2" )

