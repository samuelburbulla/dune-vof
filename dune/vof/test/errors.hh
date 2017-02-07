#ifndef DUNE_VOF_ERRORS_HH
#define DUNE_VOF_ERRORS_HH

//- C includes
#include <cmath>

//- C++ includes
#include <numeric>

//- Dune includes
#include <dune/geometry/quadraturerules.hh>

#include "interpolation.hh"


namespace Dune
{
  namespace VoF
  {

    template< class GridView, class RS, class Flags, class F >
    double l1error ( const GridView& gridView, const RS &reconstructionSet, const Flags& flags, const F &f, const double time = 0.0 )
    {
      if( gridView.comm().rank() == 0 )
        std::cout << " -- error using quadrature" << std::endl;

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
          else if ( flags.isMixed( entity ) )
            w = ( reconstructionSet[ entity ].levelSet( x ) > 0.0 );
          else
            w = 0.0;

          error += std::abs( v - w ) * qp.weight() * entity.geometry().integrationElement( qp.position() );
        }
      }

      return error;
    }


    template< class GridView, class RS, class Flags >
    double l1error ( const GridView &gridView, const RS &reconstructionSet, const Flags &flags, RotatingCircle< double, 2 > circle, const double time = 0.0 )
    {
      if( gridView.comm().rank() == 0 )
        std::cout << " -- error using intersection" << std::endl;

      double l1Error = 0.0;

      ColorFunction< GridView > uhExact( gridView );
      circleInterpolation( circle.center( time ), circle.radius( time ), uhExact );

      for ( auto entity : elements( gridView, Partitions::interior ) )
      {
        double volume = entity.geometry().volume();

        double T1 = uhExact[ entity ] * volume;
        double T0 = volume - T1;

        if ( !reconstructionSet[ entity ] )
        {
          double c = ( flags.isFull( entity ) ) ? 1.0 : 0.0;
          l1Error += T0 * std::abs( c ) + T1 * std::abs( 1.0 - c );
        }
        else
        {
          auto entityAsPolytope = makePolytope( entity.geometry() );
          auto calculatedOnePart = intersect( entityAsPolytope, reconstructionSet[ entity ], eager );
          double sharedOneVolume = intersectionVolume( calculatedOnePart, circle.center( time ), circle.radius( time ) );

          l1Error += ( T1 - sharedOneVolume ) + ( calculatedOnePart.volume() - sharedOneVolume );
        }
      }

      return l1Error;
    }


    template< class DF >
    double cellwiseL1error ( const DF& uh, RotatingCircle< double, 2 > circle, const double time = 0.0 )
    {
      if( uh.gridView().comm().rank() == 0 )
        std::cout << " -- error using cell values" << std::endl;

      DF uhExact( uh.gridView() );
      circleInterpolation( circle.center( time ), circle.radius( time ), uhExact );

      double error = 0.0;
      for ( const auto& entity : elements( uh.gridView(), Dune::Partitions::interior ) )
        error += std::abs( uh[ entity ] - uhExact[ entity ] ) * entity.geometry().volume();

      return error;
    }

    template< class CU, class F, class R, class P, class DF >
    static inline double curvatureError ( const CU &curvature, const F &flags, const R &reconstructions, const P &problem, DF &curvatureError )
    {
      double error = 0.0;
      for ( auto entity : elements( curvatureError.gridView(), Partitions::interior ) )
      {
        if ( !flags.isMixed( entity ) )
          continue;

        if ( curvature[ entity ] == 0.0 )
          continue;

        auto polygon = makePolytope( entity.geometry() );
        auto it = intersect( polygon, reconstructions[ entity ].boundary() );
        auto interface = static_cast< typename decltype( it )::Result > ( it );
        auto point = interface.centroid();

        double localError = std::abs( problem.curvature( point ) - curvature[ entity ] );
        error += localError * interface.volume();
        curvatureError[ entity ] = localError;
      }

      return error;
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_ERRORS_HH
