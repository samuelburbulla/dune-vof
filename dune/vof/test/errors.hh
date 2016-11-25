#ifndef DUNE_VOF_ERRORS_HH
#define DUNE_VOF_ERRORS_HH

//- C includes
#include <cmath>

//- C++ includes
#include <numeric>

//- Dune includes
#include <dune/geometry/quadraturerules.hh>

//- local includes
#include "utility.hh"

namespace Dune
{
  namespace VoF
  {

    template< class GridView, class RS, class Flags, class F >
    double l1error ( const GridView& gridView, const RS &reconstructionSet, const Flags& flags, const F &f, const double time = 0.0 )
    {
      using RangeType = Dune::FieldVector< double, 1 >;

      double error = 0.0;
      for ( const auto& entity : elements( gridView, Dune::Partitions::interior ) )
      {
        const auto& geo = entity.geometry();

        const auto &quad = Dune::QuadratureRules< double, GridView::dimension >::rule( geo.type(), 19 );
        for ( const auto& qp : quad )
        {
          const auto x = geo.global( qp.position() );
          RangeType v, w;
          f.evaluate( x, time, v );

          if ( flags.isFull( entity ) )
            w = 1.0;
          else if ( flags.isMixed( entity ) || flags.isFullAndMixed( entity ) )
            w = ( reconstructionSet[ entity ].levelSet( x ) > 0.0 );
          else
            w = 0.0;

          error += std::abs( v - w ) * qp.weight() * entity.geometry().integrationElement( qp.position() );
        }
      }

      return error;
    }


    template< class DF, class F, class R, class Coord >
    static inline double exactL1Error ( const DF &uhComp, const F &flags, const R &reconstructions, const Coord &center, const double radius )
    {
      double l1Error = 0.0;

      DF uhExact( uhComp );
      circleInterpolation( center, radius, uhExact );

      for ( auto entity : elements( uhComp.gridView(), Partitions::interior ) )
      {
        double volume = entity.geometry().volume();

        double T1 = uhExact[ entity ] * volume;
        double T0 = volume - T1;

        if ( !reconstructions[ entity ] )
          l1Error += T0 * std::abs( uhComp[ entity ] ) + T1 * std::abs( 1.0 - uhComp[ entity ] );
        else
        {
          auto entityAsPolytope = makePolytope( entity.geometry() );
          Dune::VoF::Polygon< Coord > calculatedOnePart = intersect( entityAsPolytope, reconstructions[ entity ] );
          double sharedOneVolume = intersectionVolume( calculatedOnePart, center, radius );

          l1Error += ( T1 - sharedOneVolume ) + ( uhComp[ entity ] * volume - sharedOneVolume );
        }
      }

      return l1Error;
    }


    template< class CU, class F, class R, class P, class DF >
    static inline double curvatureError ( const CU &curvature, const F &flags, const R &reconstructions, const P &problem, DF &curvatureError )
    {
      double error = 0.0;

      for ( auto entity : elements( curvature.gridView(), Partitions::interior ) )
      {
        if ( !flags.isMixed( entity ) && !flags.isFullAndMixed( entity ) )
          continue;

        auto polygon = makePolytope( entity.geometry() );
        auto it = intersect( polygon, reconstructions[ entity ].boundary() );
        auto interface = static_cast< typename decltype( it )::Result > ( it );

        double localError = interface.volume() * std::abs( problem.curvature( interface.centroid() ) - curvature[ entity ] );
        error += localError;
        curvatureError[ entity ] = localError;
      }

      return error;
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_ERRORS_HH
