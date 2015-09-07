
add_executable("dune-vof" main.cc polygon.hh vtu.hh)
target_link_dune_default_libraries("dune-vof")
