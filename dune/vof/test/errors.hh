#ifndef DUNE_VOF_ERRORS_HH
#define DUNE_VOF_ERRORS_HH

//- C includes
#include <cmath>

//- C++ includes
#include <numeric>

//- local includes
#include "average.hh"

namespace Dune
{
  namespace VoF
  {

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


    template< class CU, class F, class R >
    static inline double curvatureError ( const CU &curvature, const F &flags, const R &reconstructions, const double radius )
    {
      double error = 0.0;
      double kappa = 1.0 / radius;

      for ( auto entity : elements( curvature.gridView(), Partitions::interior ) )
      {
        if ( !flags.isMixed( entity ) && !flags.isFullAndMixed( entity ) )
          continue;

        auto polygon = makePolytope( entity.geometry() );
        auto it = intersect( polygon, reconstructions[ entity ].boundary() );
        auto interface = static_cast< typename decltype( it )::Result > ( it );

        error += interface.volume() * std::abs( kappa - curvature[ entity ] );
      }

      return error / ( 2.0 * M_PI * radius );
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_ERRORS_HH
