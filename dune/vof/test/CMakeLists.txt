# copy data to build source to make tests work
set( TESTDATAFILES 1dgrid.dgf 2dgrid.dgf 3dgrid.dgf parameter )
foreach( file ${TESTDATAFILES} )
  configure_file( ${file} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY )
endforeach()

dune_add_test( NAME test-geometry-2d SOURCES test-geometry.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=2" )
dune_add_test( NAME test-geometry-3d SOURCES test-geometry.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=3" )

dune_add_test( NAME test-vof-consistence SOURCES test-vof-consistence.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=2" )

dune_add_test( NAME test-curvature-2d SOURCES test-curvature.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=2" )
dune_add_test( NAME test-curvature-3d SOURCES test-curvature.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=3" )

dune_add_test( NAME test-interfacegrid-2d SOURCES test-interfacegrid.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=2" )
dune_add_test( NAME test-interfacegrid-3d SOURCES test-interfacegrid.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=3" )

dune_add_test( NAME test-vof-2d-linear SOURCES test-vof.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=2;PROBLEM=LinearWall" )
dune_add_test( NAME test-vof-2d-circle SOURCES test-vof.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=2;PROBLEM=RotatingCircle" )

dune_add_test( NAME test-vof-3d-linear SOURCES test-vof.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=3;PROBLEM=LinearWall" )
dune_add_test( NAME test-vof-3d-circle SOURCES test-vof.cc COMPILE_DEFINITIONS "SPGRID;GRIDDIM=3;PROBLEM=RotatingCircle" )
