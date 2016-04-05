
# copy data to build source to make tests work
configure_file( 3dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/3dgrid.dgf COPYONLY )
configure_file( parameter ${CMAKE_CURRENT_BINARY_DIR}/parameter COPYONLY )

add_executable( "vof-fem" main.cc ${HEADERS} )
target_link_dune_default_libraries( "vof-fem" )

set_property( TARGET vof-fem APPEND PROPERTY
  COMPILE_DEFINITIONS "ALUGRID_CUBE;GRIDDIM=3" )




add_executable( "bin2vtk" bin2vtk.cc ${HEADERS} )
target_link_dune_default_libraries( "bin2vtk" )

set_property( TARGET bin2vtk APPEND PROPERTY
COMPILE_DEFINITIONS "ALUGRID_CUBE;GRIDDIM=3" )
