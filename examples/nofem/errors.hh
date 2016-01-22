#ifndef ERRORS_HH
#define ERRORS_HH

//- C includes
#include <cmath>

//- C++ includes
#include <numeric>

//- local includes
#include "average.hh"


template< class DF, class F >
double l1error ( const DF &u, const F &f )
{
  DF solution( u.gridView() );
  average( solution, f );

  auto elements = Dune::elements( u.gridView() );
  return std::accumulate( elements.begin(), elements.end(), 0.0,
                          [ &u, &solution ] ( auto error, const auto& entity ) {
                            return error + entity.geometry().volume() * std::abs( u[ entity ] - solution [ entity ] );
                        } );
}

template< class DF, class F >
double l2error ( const DF &u, const F &ft )
{
  DF solution( u.gridView() );
  average( solution, ft );

  auto elements = Dune::elements( u.gridView() );
  return std::sqrt( std::accumulate( elements.begin(), elements.end(), 0.0,
                                     [ &u, &solution ] ( auto error, const auto &entity ) {
                                      auto val = u[ entity ] - solution[ entity ];
                                      return error + entity.geometry().volume() * val * val;
                                     }) );
}

#endif // #ifndef ERRORS_HH__
