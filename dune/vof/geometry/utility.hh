#ifndef DUNE_VOF_GEOMETRY_UTILITY_HH
#define DUNE_VOF_GEOMETRY_UTILITY_HH

#include <cassert>

#include <limits>
#include <vector>

#include <dune/common/deprecated.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/vof/brents.hh>
#include <dune/vof/geometry/halfspace.hh>
#include <dune/vof/geometry/polytope.hh>

namespace Dune
{
  namespace VoF
  {

    // distance
    // --------

    template< class ObjectA, class ObjectB >
    double distance ( const ObjectA&, const ObjectB& )
    {
      static_assert( AlwaysFalse< ObjectA >::value, "Not yet implemented!");
    };

     // hyperplane / point distance
    template< class Coord >
    double distance ( const HalfSpace< Coord >& plane, const Coord& point )
    {
      return plane.normal() * point + plane.distance();
    }

    template< class Coord >
    double distance ( const Coord& point, const HalfSpace< Coord >& plane )
    {
      return plane.normal() * point + plane.distance();
    }

    template< class DomainVector >
    inline static DomainVector rotateCCW ( const DomainVector &v )
    {
      return DomainVector{ -v[ 1 ], v[ 0 ] };
    }

     // polygon / half space   intersection (incomplete - part 1)
    template< class DomainVector, template <class> class HalfSpace >
    void polygonLineIntersection ( const Polygon2D< DomainVector > &polygon, const HalfSpace< DomainVector > &g, Polygon2D< DomainVector > &intersectionPolygon )
    {
      using std::abs;

      for ( std::size_t i = 0; i < polygon.corners(); ++i )
      {
        auto ci  = distance( g, polygon[ i ] );
        auto cip = distance( g, polygon[ i+1 ] );

        if( abs( ci ) < 1e-12 )
          intersectionPolygon.addVertex( polygon[i] );
        else if( ( ci > 0.0 ) ^ ( cip > 0.0 ) )
        {
          DomainVector point;
          point.axpy( -cip / ( ci - cip ), polygon[ i ] );
          point.axpy(  ci  / ( ci - cip ), polygon[ i+1 ] );

          // add intersection point
          intersectionPolygon.addVertex( point );
        }
      }
    }

     // cell / line intersection
    template< class Geo, template <class> class HalfSpace, class DomainVector >
    void lineCellIntersections ( const Geo &geo,
                                 const HalfSpace< DomainVector > &g,
                                 std::vector< DomainVector > &intersections,
                                 const double TOL = 1e-12
      )
    {
      using std::abs;

      const int dim = 2;  // only two-dimensional
      intersections.clear();
      const auto &refElement = Dune::ReferenceElements< double, dim >::general( geo.type() );

      for( int k = 0; k < refElement.size( dim-1 ); ++k )
      {
        int i = refElement.subEntity( k, dim-1, 0, dim );
        int j = refElement.subEntity( k, dim-1, 1, dim );

        const DomainVector& x0 = geo.corner( i );
        const DomainVector& x1 = geo.corner( j );

        auto c0 = distance( g, x0 );
        auto c1 = distance( g, x1 );

        if( ( c0 > 0.0 ) ^ ( c1 > 0.0 ) )
        {
          DomainVector point;
          point.axpy( -c1 / ( c0 - c1 ), x0 );
          point.axpy(  c0 / ( c0 - c1 ), x1 );

          // add intersection point
          intersections.push_back( point );
        }
        else if( abs( c0 ) < TOL  )
          intersections.push_back( x0 );
        else if( abs( c1 ) < TOL )
          intersections.push_back( x1 );
      }

      std::sort( intersections.begin(), intersections.end(),
        [] ( const DomainVector& v, const DomainVector& w )
        {
          if ( v[ 0 ] == w[ 0 ] )
            return v[ 1 ] < w[ 1 ];
          return v[ 0 ] < w[ 0 ];
        } );

      auto it = std::unique( intersections.begin(), intersections.end(),
        [] ( const DomainVector& v, const DomainVector& w )
        {
          return ( v - w ).two_norm2() < std::numeric_limits< typename DomainVector::value_type >::epsilon();
        } );
      intersections.resize( std::distance( intersections.begin(), it ) );
    }

     // polygon / half space   intersection (incomplete - part 2)
    template< class DomainVector, class HalfSpace >
    void polyAddInnerVertices ( const Polygon2D< DomainVector > &sourcePolygon, const HalfSpace &g, Polygon2D< DomainVector >& endPolygon, const double TOL = 1e-12 )
    {
      for( std::size_t i = 0; i < sourcePolygon.corners(); ++i )
        if( distance( sourcePolygon[ i ], g ) > 0.0 )
          endPolygon.addVertex( sourcePolygon[ i ] );
    }


    template< class Geo, template <class> class HalfSpace, class DomainVector >
    double getVolumeFraction ( const Geo &geo, const HalfSpace< DomainVector > &g )
    {
      Polygon2D< DomainVector > polygon;

      for( int i = 0; i < geo.corners(); ++i )
        if( distance( geo.corner( i ), g ) > 0.0 )
          polygon.addVertex( geo.corner( i ) );

      std::vector< DomainVector > intersections;
      lineCellIntersections( geo, g, intersections );
      for( auto v : intersections )
        polygon.addVertex( v );

      return polygon.volume() / geo.volume();
    }


    template< class Geometry, class Color, class Reconstruction, class PointList>
    void computeInterfaceLinePosition ( const Geometry &geo, const Color &concentration, Reconstruction &g, PointList &intersections )
    {
      using vtype = decltype( geo.volume() );
      using limits = std::numeric_limits< vtype >;
      using ctype = typename Reconstruction::ctype;

      ctype pMin = 0, pMax = 0;
      vtype volume,
            volMin = limits::lowest(), // lower bound
            volMax = limits::max(); // upper bound

      // upper/lower bound for
      for( int i = 0; i < geo.corners(); ++i )
      {
        Reconstruction h( g.normal(), geo.corner( i ) );

        volume = getVolumeFraction( geo, h );

        if( ( volume <= volMax ) && ( volume >= concentration ) )
        {
          pMax = h.distance();
          volMax = volume;
        }
        if( ( volume >= volMin ) && ( volume <= concentration ) )
        {
          pMin = h.distance();
          volMin = volume;
        }
      }
      assert( volMin <= volMax );

      if ( volMax == limits::max() )
        g.distance() = pMin;
      else if ( volMin == limits::lowest() )
        g.distance() = pMax;
      else if ( volMin == volMax )
        g.distance() = pMin;
      else
        g.distance() = brentsMethod( [ &geo, &concentration, &g ] ( ctype p ) -> ctype {
                                        Reconstruction h( g.normal(), p );
                                        return ( getVolumeFraction( geo, h ) - concentration );
                                     }, pMin, pMax, 1e-12 );

      lineCellIntersections( geo, g, intersections );
    }


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_UTILITY_HH
