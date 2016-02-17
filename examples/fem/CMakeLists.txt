
# copy data to build source to make tests work
configure_file( 2dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY )
configure_file( parameter ${CMAKE_CURRENT_BINARY_DIR}/parameter COPYONLY )
configure_file( allbin2vtk.sh ${CMAKE_CURRENT_BINARY_DIR}/allbin2vtk.sh COPYONLY )

add_executable( "vof-fem" main.cc ${HEADERS} )
target_link_dune_default_libraries( "vof-fem" )

set_property( TARGET vof-fem APPEND PROPERTY
  COMPILE_DEFINITIONS "ALUGRID_CUBE;GRIDDIM=2" )




add_executable( "bin2vtk" bin2vtk.cc ${HEADERS} )
target_link_dune_default_libraries( "bin2vtk" )

set_property( TARGET bin2vtk APPEND PROPERTY
COMPILE_DEFINITIONS "ALUGRID_CUBE;GRIDDIM=2" )
