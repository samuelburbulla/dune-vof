#ifndef DUNE_VOF_GEOMETRICUTILITY_HH
#define DUNE_VOF_GEOMETRICUTILITY_HH

#include <cassert>

#include <limits>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/vof/brents.hh>

namespace Dune
{
  namespace VoF
  {

    template< class DomainVector >
    inline static DomainVector rotateCCW ( const DomainVector &v )
    {
      return DomainVector{ -v[ 1 ], v[ 0 ] };
    }

    template< template <class> class Hyperplane, class DomainVector >
    const DomainVector lineIntersection ( const Hyperplane< DomainVector > &g, const Hyperplane< DomainVector > &l )  // make dim-universial
    {

      assert ( DomainVector::dimension == 2 );

      Dune::FieldMatrix< double, 2, 2 > A( { g.normal(), l.normal() }  );
      DomainVector b( { -g.distance(), -l.distance() } );
      DomainVector x;

      A.solve( x, b );

      return x;
    }


    // konvexes 2D-Polygon
    template< class DomainVector >
    struct Polygon2D
    {
      Polygon2D< DomainVector >() {}

      const DomainVector operator[] ( const int i ) const { return points[ i % points.size() ]; }

      void addVertex ( const DomainVector &vertex, const bool correctOrder = false, const double TOL = 1e-12 )
      {
        std::size_t n = points.size();

        // no specific order needed or vertex explicitly given in correct order
        if( n == 0 || n == 1 || correctOrder )
        {
          points.push_back( vertex );
        }

        // insert new vertex in counterclockwise order
        for( std::size_t i = 0; i < n; ++i )
        {
          auto normal = points[ (i+1)%n ];
          normal -= points[ i ];
          rotccw( normal );

          auto center = points[ i ];
          center += points[ (i+1)%n ];
          center *= 0.5;

          center -= vertex;

          if( normal * center > 0 )
          {
            points.insert( points.begin() + i + 1, vertex );
            return;
          }
        }
      }

      const std::size_t corners () const { return points.size(); }

      double volume () const
      {
        double sum = 0;
        int n = points.size();
        for( int i = 0; i < n; i++ )
          sum += ( points[ i ][ 1 ] + points[ (i+1)%n ][ 1 ] ) * ( points[ i ][ 0 ] - points[ (i+1)%n ][ 0 ] );
        return sum / 2.0;
      }

      bool pointInBorders ( const DomainVector &vertex, const double TOL = 1e-12 ) const
      {
        int n = corners();

        if( n == 0 )
          return false;

        for( int i = 0; i < n; ++i )
        {
          auto edge = points[ (i+1)%n ];
          edge -= points[ i ];
          rotccw ( edge );

          auto skalar = edge * ( vertex - points[ i ] );
          if (  skalar < 0 && std::abs( skalar ) > TOL ) return false;
        }

        return true;
      }

      void clear () { points.clear(); }

    private:
      std::vector< DomainVector > points;
    };

    template< class DomainVector, template <class> class Hyperplane >
    void polygonLineIntersection ( const Polygon2D< DomainVector > &polygon, const Hyperplane< DomainVector > &g, Polygon2D< DomainVector > &intersectionPolygon )
    {
      for ( std::size_t i = 0; i < polygon.corners(); ++i )
      {
        if( isOnRecLine( polygon[i], g ) )
        {
          intersectionPolygon.addVertex( polygon[i] );
        }
        else if( isInner( polygon[i], g ) ^ isInner( polygon[i+1], g ) )
        {
          auto normal = polygon[ i ];
          normal -= polygon[ i+1 ];

          const Hyperplane< DomainVector > lineThroughEdge( rotateCCW( normal ), polygon[ i ] );

          // add intersection point
          intersectionPolygon.addVertex( lineIntersection( g, lineThroughEdge ) );
        }
      }
    }

    template< class DomainVector >
    bool dvComp ( const DomainVector &v, const DomainVector &w ) {
      if ( v[0] == w[0] ) return v[1] < w[1];
      return v[0] < w[0];
    }

    template< class DomainVector >
    bool dvEq ( const DomainVector &v, const DomainVector &w ) {
      return ( ( v - w ).one_norm() < 1e-12 );
    }

    template< class Geo, template <class> class Hyperplane, class DomainVector >
    void lineCellIntersections (
      const Geo &geo,
      const Hyperplane< DomainVector > &g,
      std::vector< DomainVector > &intersections,
      const double TOL = 1e-12
      )
    {
      const int dim = 2;  // only two-dimensional
      intersections.clear();

      const auto &refElement = Dune::ReferenceElements< double, dim >::general( geo.type() );

      for( int k = 0; k < refElement.size( dim-1 ); ++k )
      {
        int i = refElement.subEntity( k, dim-1, 0, dim );
        int j = refElement.subEntity( k, dim-1, 1, dim );

        const DomainVector& c0 = geo.global( refElement.position( i, dim ) );
        const DomainVector& c1 = geo.global( refElement.position( j, dim ) );

        if( isInner( c0, g, TOL ) ^ isInner( c1, g, TOL ) )
        {
          // build line through edge for intersection
          auto normal = c0;
          normal -= c1;
          const Hyperplane< DomainVector > lineThroughEdge( rotateCCW( normal ), c0 );

          // add intersection point
          intersections.push_back( lineIntersection( g, lineThroughEdge ) );
        }
        else if( isOnRecLine( c0, g, TOL ) ) intersections.push_back( c0 );
        else if( isOnRecLine( c1, g, TOL ) ) intersections.push_back( c1 );
      }

      std::sort( intersections.begin(), intersections.end(), dvComp< DomainVector > );
      auto it = std::unique( intersections.begin(), intersections.end(), dvEq< DomainVector > );
      intersections.resize( std::distance( intersections.begin(), it ) );
    }


    template < template <class> class Hyperplane, class DomainVector >
    bool isInner ( const DomainVector &vertex, const Hyperplane< DomainVector > &g, const double TOL = 1e-12 )
    {
      return ( vertex * g.normal() + g.distance() ) >= TOL;
    }

    template < template <class> class Hyperplane, class DomainVector >
    bool isOnRecLine ( const DomainVector &vertex, const Hyperplane< DomainVector > &g, const double TOL = 1e-12 )
    {
      return std::abs( vertex * g.normal() + g.distance() ) < TOL;
    }


    template< class DomainVector, class Hyperplane >
    void polyAddInnerVertices ( const Polygon2D< DomainVector > &sourcePolygon, const Hyperplane &g, Polygon2D< DomainVector >& endPolygon, const double TOL = 1e-12 )
    {
      for( std::size_t i = 0; i < sourcePolygon.corners(); ++i )
        if( isInner( sourcePolygon[ i ], g ) )
          endPolygon.addVertex( sourcePolygon[ i ] );
    }

    template< class Geo, template <class> class Hyperplane, class DomainVector >
    double getVolumeFraction ( const Geo &geo, const Hyperplane< DomainVector > &g )
    {
      Polygon2D< DomainVector > polygon;

      for( int i = 0; i < geo.corners(); ++i )
        if( isInner( geo.corner( i ), g ) )
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

#endif // #ifndef DUNE_VOF_GEOMETRICUTILITY_HH
