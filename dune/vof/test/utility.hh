#ifndef DUNE_VOF_UTILITY_HH
#define DUNE_VOF_UTILITY_HH

// dune-fem includes
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/gridpart/leafgridpart.hh>

// dune-vof includes
#include <dune/vof/geometry/2d/polygon.hh>
#include <dune/vof/geometry/intersect.hh>

namespace Dune
{
  namespace VoF
  {
    // Exact interpolation for circle
    // ==============================


    template < class Coordinate >
    static inline int circleIntersection ( const Dune::VoF::Line< Coordinate >& line, const Coordinate& c, double r, std::vector< Coordinate > &points )
    {
      points.clear();

      // Quelle: Wikipedia ( https://de.wikipedia.org/wiki/Schnittpunkt )
      Coordinate n = generalizedCrossProduct( line.vertex(1) - line.vertex(0) );
      double p = n * line.vertex( 0 );
      double d = p - n[0] * c[0] - n[1] * c[1];

      double diskr = r * r * n.two_norm2() - d * d;
      if ( diskr < 0 )
        return 0;
      else if ( std::abs( diskr ) == 0.0 )
      {
        Coordinate v ( { n[0] * d, n[1] * d } );
        v /= n.two_norm2();
        v += c;
        if ( ( line.vertex(0) - v ) * ( line.vertex(1) - v ) <= 0.0 )
        {
          points.push_back( v );
          return 1;
        }
        return 0;
      }
      else
      {
        double sqrtDiskr = std::sqrt( diskr );
        Coordinate v1 ( { n[0] * d + n[1] * sqrtDiskr, n[1] * d - n[0] * sqrtDiskr  } );
        v1 /= n.two_norm2();
        v1 += c;
        if ( ( line.vertex(0) - v1 ) * ( line.vertex(1) - v1 ) <= 0.0 )
          points.push_back( v1 );

        Coordinate v2 ( { n[0] * d - n[1] * sqrtDiskr, n[1] * d + n[0] * sqrtDiskr  } );
        v2 /= n.two_norm2();
        v2 += c;
        if ( ( line.vertex(0) - v2 ) * ( line.vertex(1) - v2 ) <= 0.0 )
          points.push_back( v2 );

        // Check order of points
        if ( points.size() == 2 )
        {
          Dune::VoF::Line< Coordinate > lineBetweenPoints ( points[0], points[1] );
          Coordinate n1 = generalizedCrossProduct( lineBetweenPoints.vertex(1) - lineBetweenPoints.vertex(0) );
          if ( n1 * n < 0 )
            std::swap( points[0], points[1] );
        }

        return 2;
      }

    }


    template< class Polygon, class Coordinate >
    static inline double intersectionVolume ( const Polygon &polygon, const Coordinate &center, double radius )
    {
      std::vector< Coordinate > secPoints;
      std::vector< Coordinate > polyPoints;

      // Start with vertex outside of circle
      int i = 0, i0 = 0;
      for (; i < polygon.size(); ++i )
        if ( ( polygon.vertex( i ) - center ).two_norm() > radius )
        {
          i0 = i;
          break;
        }
      // otherwise, entity is completely included in circle
      if ( i == polygon.size() )
        return polygon.volume();

      for ( int i = 0; i < polygon.size(); ++i )
      {
        auto edge = polygon.edge( ( i + i0 ) % polygon.size() );
        std::vector< Coordinate > intersections;
        circleIntersection ( edge, center, radius, intersections );

        switch ( intersections.size() )
        {
          case 1:
            secPoints.push_back( intersections[ 0 ] );
            polyPoints.push_back( intersections[ 0 ] );
            break;
          case 2:
            secPoints.push_back( intersections[ 0 ] );
            secPoints.push_back( intersections[ 1 ] );
            polyPoints.push_back( intersections[ 0 ] );
            polyPoints.push_back( intersections[ 1 ] );
            break;
          default:
            break;
        };

        if ( ( edge.vertex( 1 ) - center ).two_norm() < radius )
          polyPoints.push_back( edge.vertex( 1 ) );
      }

      double volume = 0.0;

      if ( secPoints.size() > 1 )
      {
        if ( polyPoints.size() > 2 )
        {
          Dune::VoF::Polygon< Coordinate > poly ( polyPoints );
          volume += poly.volume();
        }

        for ( std::size_t i = 1; i < secPoints.size(); i = i+2 )
        {
          // Quelle: Wikipedia (Kreissegment)
          Dune::VoF::Line< Coordinate > sec ( secPoints[i], secPoints[ (i+1) % secPoints.size() ] );
          double alpha = 2.0 * std::asin( sec.volume() / ( 2.0 * radius ) );
          volume += radius * radius / 2.0 * ( alpha - std::sin( alpha ) );
        }
      }

      return volume;
    }


    template< class Coord, class DF >
    static inline void circleInterpolation ( const Coord &center, double radius, DF &uh )
    {
      uh.clear();

      for ( const auto& entity : elements( uh.gridView(), Partitions::interior ) )
      {
        const auto& geo = entity.geometry();
        Dune::VoF::Polygon< Coord > polygon = Dune::VoF::makePolytope( geo );

        uh[ entity ] = intersectionVolume( polygon, center, radius ) / geo.volume();
      }
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



  } // namespace VoF

}  // namespace Dune

#endif

