set(HEADERS
  cartesianheightfunctioncurvature.hh
  circleinterpolationcurvature.hh
  generalheightfunctioncurvature.hh
  interfacenormaldivergencecurvature.hh
  normaldivergencecurvature.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof/curvature)
