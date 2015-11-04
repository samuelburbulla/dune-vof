configure_file( parameter.ini ${CMAKE_CURRENT_BINARY_DIR}/parameter.ini COPYONLY )
configure_file( cube2.dgf ${CMAKE_CURRENT_BINARY_DIR}/cube2.dgf COPYONLY )

add_executable(
  "dune-vof" 
  main.cc 
  domain.hh
  testproblem.hh
  colorFunction.hh
)
target_link_dune_default_libraries("dune-vof")
set_property( TARGET dune-vof APPEND PROPERTY
  COMPILE_DEFINITIONS "ALUGRID_SIMPLEX;GRIDDIM=2" )
