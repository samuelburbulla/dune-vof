#install headers
set(vof_HEADERS
  domain.hh
  errors.hh
  evolve.hh
  flagCells.hh
  geometricaltoolbox.hh
  initialize.hh
  reconstruct.hh
  testproblem.hh
  vof.hh
  vtkout.hh)
install(FILES ${vof_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof)
