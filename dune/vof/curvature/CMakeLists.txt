set( curvaturedir ${CMAKE_INSTALL_INCLUDEDIR}/dune/vof/curvature )
set( curvature_HEADERS
  cartesianheightfunctioncurvature.hh
  circleinterpolationcurvature.hh
  generalheightfunctioncurvature.hh
  interfacenormaldivergencecurvature.hh
  normaldivergencecurvature.hh
)

install(FILES ${curvature_HEADERS} DESTINATION ${curvaturedir})
