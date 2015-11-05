#ifndef ERRORS_HH
#define ERRORS_HH

//- dune-geometry includes
#include <dune/geometry/quadraturerules.hh>

//- local includes
#include "average.hh"


template< class DF, class F >
double l1error ( const DF &u, const F &f )
{
  DF solution( u.gridView() );
  average( solution, f );

  double error = 0;
  for( const auto& entity : Dune::elements( u.gridView() ) )
  {
    const auto geo = entity.geometry();

    error += geo.volume() * std::abs( u[ entity ] - solution[ entity ] );
  }

  return error;
}

template< class DF, class F >
double l2error ( const DF &u, const F &ft )
{
  DF solution( u.gridView() );
  average( solution, ft );

  double error = 0;

  for( const auto& entity : Dune::elements( u.gridView() ) )
  {
    const auto geo = entity.geometry();

    double diff = u[ entity ] - solution[ entity ];
    diff *= diff;
    error += geo.volume() * diff;
  }

  return std::sqrt( error );
}

#endif // #ifndef ERRORS_HH__
