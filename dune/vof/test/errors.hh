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

#endif // #ifndef ERRORS_HH__
