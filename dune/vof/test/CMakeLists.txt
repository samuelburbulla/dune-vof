# if no grid type is specified during configuration time we choose yaspgrid
if( NOT GRIDTYPE )
  set( GRIDTYPE YASPGRID )
endif()

# same as for GridType, we set the default grid dimension to 2
if( NOT GRIDDIM )
  set( GRIDDIM 2 )
endif()

# copy data to build source to make tests work
set( TESTDATAFILES 2dgrid.dgf 3dgrid.dgf parameter.ini )
foreach( file ${TESTDATAFILES} )
  configure_file( ${file} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY )
endforeach()

# some default flags needed very often
set( DEFAULTFLAGS "COUNT_FLOPS;${GRIDTYPE};GRIDDIM=${GRIDDIM}" )

# variable which gathers all test targets
set( TESTS )

dune_add_test( NAME test-vof-2d-circle SOURCES test-vof-2d-circle.cc COMPILE_DEFINITIONS "${DEFAULTFLAGS}" )
