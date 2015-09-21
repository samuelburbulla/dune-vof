configure_file( parameter.ini ${CMAKE_CURRENT_BINARY_DIR}/parameter.ini COPYONLY )

add_executable("dune-vof" main.cc polygon.hh vtu.hh)
target_link_dune_default_libraries("dune-vof")
