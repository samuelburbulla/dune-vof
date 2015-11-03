configure_file( parameter.ini ${CMAKE_CURRENT_BINARY_DIR}/parameter.ini COPYONLY )
configure_file( cube2.dgf ${CMAKE_CURRENT_BINARY_DIR}/cube2.dgf COPYONLY )

add_executable(
  "dune-vof" 
  main.cc 
  polygon.hh 
  vtu.hh
  domain.hh
  errors.hh
  evolve.hh
  geometricaltoolbox.hh
  initialize.hh
  reconstruct.hh
  testproblem.hh
  io.hh
  secondyoungsnormalguessing.hh
  brent.hh
  cputime.hh
)
target_link_dune_default_libraries("dune-vof")
set_property( TARGET dune-vof APPEND PROPERTY
  COMPILE_DEFINITIONS "ALUGRID_SIMPLEX;GRIDDIM=2" )
