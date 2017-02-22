
set( GRIDTYPE SPGRID )

set( GRIDDIM 2 )

set( GRIDFILES
  parameter
)



# copy data to build source to make tests work
foreach( file ${GRIDFILES} )
  configure_file( ${file} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY )
endforeach()


set( EXECUTABLES
  vof
  bin2vtk
)

foreach ( exec ${EXECUTABLES} )
  add_executable( "${exec}" ${exec}.cc ${HEADERS} )
  target_link_dune_default_libraries( "${exec}" )
  set_property( TARGET ${exec} APPEND PROPERTY COMPILE_DEFINITIONS "${GRIDTYPE}" "GRIDDIM=${GRIDDIM}" LINK_LIBRARIES dunefem )
endforeach()

