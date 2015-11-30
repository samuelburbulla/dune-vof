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
    inline static void rotccw ( DomainVector &v )
    {
      auto t = v[0];
      v[0] = -v[1];
      v[1] = t;
    }


    template< template <class> class HyperSurface, class DomainVector >
    const DomainVector lineIntersection ( const HyperSurface< DomainVector > &g, const HyperSurface< DomainVector > &l )  // make dim-universial
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

      void addVertex ( const DomainVector &vertex, const double TOL = 1e-12 )
      {
        std::size_t n = points.size();

        // check, if point is already here
        for( std::size_t i = 0; i < n; i++ )
          if( (vertex - points[ i ]).two_norm() < TOL )
            return;

        // insert new vertex in counterclockwise order
        if( n == 0 || n == 1 )
        {
          points.push_back( vertex );
          return;
        }
        else
          // start with last inserted edge, for the case that points are inserted in correct order
          for( std::size_t i = n - 1; i >= 0; --i )
          {
            auto normal = points[ (i+1)%n ];
            normal -= points[ i ];
            rotccw( normal );

            auto center = points[ i ];
            center += points[ (i+1)%n ];
            center *= 0.5;

            center -= vertex;
            center *= -1.0;

            if( normal * center < 0 )
            {
              points.insert( points.begin() + i + 1, vertex );
              //std::cerr << n-i << "/" << n << std::endl;
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

      bool pointInBorders ( const DomainVector &vertex ) const
      {
        int n = points.size();

        if( n == 0 )
          return false;

        bool inside = true;

        for( int i = 0; i < n; ++i )
          inside = inside && this->SkalarProdTest( vertex, points[ i ], points[ (i+1)%n ] );

        return inside;
      }

      void clear () { points.clear(); }

      std::vector< DomainVector > points;

    private:
      const int SkalarProdTest ( const DomainVector &vertex, const DomainVector &p1, const DomainVector &p2, const double TOL = 1e-12 ) const
      {
        auto p = p2 - p1;
        rotccw ( p );

        auto skalar = ( p * ( vertex - p1) );

        return ( skalar >= 0 || std::abs( skalar ) < TOL );
      }

    };

    template< class DomainVector, template <class> class HyperSurface >
    void polygonLineIntersection ( const Polygon2D< DomainVector > &polygon, const HyperSurface< DomainVector > &g, Polygon2D< DomainVector > &intersectionPolygon )
    {
      for ( std::size_t i = 0; i < polygon.corners(); ++i )
      {
        if( isOnRecLine( polygon[i], g ) )
        {
          intersectionPolygon.addVertex( polygon[i] );
        }
        else if( isInner( polygon[i], g ) ^ isInner( polygon[i+1], g ) )
        {
          auto normal = polygon[ i ] - polygon[ i+1 ];
          rotccw ( normal );
          const HyperSurface< DomainVector > lineThroughEdge( normal, polygon[ i ] );

          // add intersection point
          intersectionPolygon.addVertex( lineIntersection( g, lineThroughEdge ) );
        }
      }
    }

    template< class DomainVector >
    void insertElementIfNotExists ( const DomainVector &v, std::vector< DomainVector > &list, const double TOL = 1e-12 )
    {

      for( std::size_t i = 0; i < list.size(); ++i )
        if( ( list[ i ] - v).two_norm() < TOL )
          return;

      list.push_back( v );
    }

    template< class Geo, template <class> class HyperSurface, class DomainVector >
    std::vector< DomainVector > lineCellIntersections ( const Geo &geo, const HyperSurface< DomainVector > &g, const double TOL = 1e-12 )
    {
      const int dim = 2;

      std::vector< DomainVector > intersectionPoints;

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
          auto normal = c0 - c1;
          rotccw( normal );
          const HyperSurface< DomainVector > lineThroughEdge( normal, c0 );

          // add intersection point
          intersectionPoints.push_back( lineIntersection( g, lineThroughEdge ) );
        }
        else if( isOnRecLine( c0, g, TOL ) )
          insertElementIfNotExists( c0, intersectionPoints );
        else if( isOnRecLine( c1, g, TOL ) )
          insertElementIfNotExists( c1, intersectionPoints );
      }

      return intersectionPoints;
    }

    template < template <class> class HyperSurface, class DomainVector >
    bool isInner ( const DomainVector &vertex, const HyperSurface< DomainVector > &g, const double TOL = 1e-12 )
    {
      return ( vertex * g.normal() + g.distance() ) >= TOL;
    }

    template < template <class> class HyperSurface, class DomainVector >
    bool isOnRecLine ( const DomainVector &vertex, const HyperSurface< DomainVector > &g, const double TOL = 1e-12 )
    {
      return std::abs( vertex * g.normal() + g.distance() ) < TOL;
    }

    template< class Geo, class DomainVector, class HyperSurface >
    void polyAddInnerVertices ( const Geo &geo, const HyperSurface &g, Polygon2D< DomainVector >& polygon, const double TOL = 1e-12 )
    {
      for( int i = 0; i < geo.corners(); ++i )
        if( isInner( geo.corner( i ), g, TOL ) )
          polygon.addVertex( geo.corner( i ) );
    }

    template< class DomainVector, class HyperSurface >
    void polyAddInnerVertices ( const Polygon2D< DomainVector >& sourcePolygon, const HyperSurface &g, Polygon2D< DomainVector >& endPolygon, const double TOL = 1e-12 )
    {
      for( std::size_t i = 0; i < sourcePolygon.corners(); ++i )
        if( isInner( sourcePolygon[ i ], g, TOL ) )
          endPolygon.addVertex( sourcePolygon[ i ] );
    }

    template< class Geo, template <class> class HyperSurface, class DomainVector >
    double getVolumeFraction ( const Geo &geo, const HyperSurface< DomainVector > &g )
    {
      Polygon2D< DomainVector > polygonVertices;

      auto lip = lineCellIntersections( geo, g );
      for( auto &v : lip )
        polygonVertices.addVertex( v );

      polyAddInnerVertices( geo, g, polygonVertices );

      assert( polygonVertices.corners() <= 5 );

      return polygonVertices.volume() / geo.volume();
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

      intersections.clear();
      intersections = lineCellIntersections( geo, g );
    }


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRICUTILITY_HH
