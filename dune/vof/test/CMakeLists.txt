# copy data to build source to make tests work
set( TESTDATAFILES 2dgrid.dgf 3dgrid.dgf parameter.ini )
foreach( file ${TESTDATAFILES} )
  configure_file( ${file} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY )
endforeach()

dune_add_test( NAME test-geometry-2d SOURCES test-geometry.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=2" )
dune_add_test( NAME test-geometry-3d SOURCES test-geometry.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=3" )

dune_add_test( NAME test-vof-2d-circle-exact SOURCES test-vof-2d-circle-exact.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=2" )

dune_add_test( NAME test-vof-3d-circle SOURCES test-vof.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=3;PROBLEM=RotatingCircle" )
dune_add_test( NAME test-vof-2d-linear SOURCES test-vof.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=2;PROBLEM=LinearWall" )
dune_add_test( NAME test-vof-3d-linear SOURCES test-vof.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=3;PROBLEM=LinearWall" )
dune_add_test( NAME test-vof-2d-diagonal SOURCES test-vof.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=2;PROBLEM=Diagonal" )
dune_add_test( NAME test-vof-3d-diagonal SOURCES test-vof.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=3;PROBLEM=Diagonal" )
