# copy data to build source to make tests work
set( TESTDATAFILES 2dgrid.dgf 3dgrid.dgf parameter.ini )
foreach( file ${TESTDATAFILES} )
  configure_file( ${file} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY )
endforeach()

dune_add_test( NAME test-vof-2d-circle SOURCES test-vof-2d-circle.cc COMPILE_DEFINITIONS "YASPGRID;GRIDDIM=2" )
dune_add_test( NAME test-vof-3d-circle SOURCES test-vof-3d-circle.cc COMPILE_DEFINITIONS "YASPGRID;GRIDDIM=3" )
dune_add_test( NAME test-vof-2d-linear SOURCES test-vof-2d-linear.cc COMPILE_DEFINITIONS "YASPGRID;GRIDDIM=2" )
# dune_add_test( NAME test-vof-3d-linear SOURCES test-vof-3d-linear.cc COMPILE_DEFINITIONS "YASPGRID;GRIDDIM=3" )
