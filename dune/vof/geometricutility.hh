#ifndef DUNE_VOF_GEOMETRICUTILITY_HH
#define DUNE_VOF_GEOMETRICUTILITY_HH

#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/vof/brents.hh>

namespace Dune
{
  namespace VoF
  {

    template< class V >
    inline static V rotate90degreesCounterClockwise ( const V &v )
    {
      return V{ -v[ 1 ], v[ 0 ] };
    }


    template< template <class> class HyperSurface, class fvector >
    const fvector lineIntersection ( const HyperSurface< fvector > &g, const HyperSurface< fvector > &l )  // make dim-universial
    {

      assert ( fvector::dimension == 2 );

      Dune::FieldMatrix< double, 2, 2 > A( { g.normal(), l.normal() }  );
      fvector b( { -g.distance(), -l.distance() } );
      fvector x;

      A.solve( x, b );

      return x;
    }


    // konvexes 2D-Polygon
    template< class V >
    struct Polygon2D
    {
      Polygon2D< V >() {}

      const V operator[] ( const int i ) const { return points[ i % points.size() ]; }

      void addVertex ( const V &vertex, const double TOL = 1e-12 )
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
          for( std::size_t i = 0; i < points.size(); i++ )
          {
            V normal = points[ (i+1)%n ];
	          normal -= points[ i ];

            if( ( rotate90degreesCounterClockwise( normal ) * ( vertex - points[ i ] ) ) < 0 )
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

      const bool pointInBorders ( const V &vertex ) const
      {
        int n = points.size();

        if( n == 0 )
          return false;

        bool inside = true;

        for( int i = 0; i < n; ++i )
          inside = inside && this->SkalarProdTest( vertex, points[ i ], points[ (i+1)%n ] );

        return inside;
      }

      std::vector< V > points;
      double vol = 0;

    private:
      const int SkalarProdTest ( const V &vertex, const V &p1, const V &p2, const double TOL = 1e-12 ) const
      {
        auto skalar = ( rotate90degreesCounterClockwise( p2 - p1 ) * ( vertex - p1) );

        return ( skalar >= 0 || std::abs( skalar ) < TOL );
      }

    };

    template< class V, template <class> class HyperSurface >
    void polygonLineIntersection ( const Polygon2D< V > &polygon, const HyperSurface< V > &g, Polygon2D< V > &intersectionPolygon )
    {
      for ( std::size_t i = 0; i < polygon.corners(); ++i )
      {
        if( isOnRecLine( polygon[i], g ) )
        {
          intersectionPolygon.addVertex( polygon[i] );
        }
        else if( isInner( polygon[i], g ) ^ isInner( polygon[i+1], g ) )
        {
          const HyperSurface< V > lineThroughEdge( rotate90degreesCounterClockwise( polygon[ i ] - polygon[ i+1 ] ), polygon[ i ] );

          // add intersection point
          intersectionPolygon.addVertex( lineIntersection( g, lineThroughEdge ) );
        }
      }
    }

    template< class V >
    void insertElementIfNotExists ( const V &v, std::vector< V > &list, const double TOL = 1e-12 )
    {

      for( std::size_t i = 0; i < list.size(); ++i )
        if( ( list[ i ] - v).two_norm() < TOL )
          return;

      list.push_back( v );
    }

    template< class Geo, template <class> class HyperSurface, class V >
    std::vector< V > lineCellIntersections ( const Geo &geo, const HyperSurface<V> &g, const double TOL = 1e-12 )
    {
      const int dim = 2;

      std::vector< V > intersectionPoints;

      const auto &refElement = Dune::ReferenceElements< double, dim >::general( geo.type() );

      for( int k = 0; k < refElement.size( dim-1 ); ++k )
      {
        int i = refElement.subEntity( k, dim-1, 0, dim );
        int j = refElement.subEntity( k, dim-1, 1, dim );

        const V& c0 = geo.global( refElement.position( i, dim ) );
        const V& c1 = geo.global( refElement.position( j, dim ) );

        if( isInner( c0, g, TOL ) ^ isInner( c1, g, TOL ) )
        {
          // build line through edge for intersection
          const HyperSurface< V > lineThroughEdge( rotate90degreesCounterClockwise( c0 - c1 ), c0 );

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

    template < template <class> class HyperSurface, class V >
    bool isInner ( const V &vertex, const HyperSurface<V> &g, const double TOL = 1e-12 )
    {
      return ( vertex * g.normal() + g.distance() ) >= TOL;
    }

    template < template <class> class HyperSurface, class V >
    bool isOnRecLine ( const V &vertex, const HyperSurface<V> &g, const double TOL = 1e-12 )
    {
      return std::abs( vertex * g.normal() + g.distance() ) < TOL;
    }

    template< class Geo, class V, class HyperSurface >
    void polyAddInnerVertices ( const Geo &geo, const HyperSurface &g, Polygon2D< V >& polygon, const double TOL = 1e-12 )
    {
      for( int i = 0; i < geo.corners(); ++i )
        if( isInner( geo.corner( i ), g, TOL ) )
          polygon.addVertex( geo.corner( i ) );
    }

    template< class V, class HyperSurface >
    void polyAddInnerVertices ( const Polygon2D< V >& sourcePolygon, const HyperSurface &g, Polygon2D< V >& endPolygon, const double TOL = 1e-12 )
    {
      for( std::size_t i = 0; i < sourcePolygon.corners(); ++i )
        if( isInner( sourcePolygon[ i ], g, TOL ) )
          endPolygon.addVertex( sourcePolygon[ i ] );
    }

    template< class Geo, template <class> class HyperSurface, class V >
    double getVolumeFraction ( const Geo &geo, const HyperSurface<V> &g )
    {
      Polygon2D< V > polygonVertices;

      auto lip = lineCellIntersections( geo, g );
      for( auto &v : lip )
        polygonVertices.addVertex( v );

      polyAddInnerVertices( geo, g, polygonVertices );

      assert( polygonVertices.corners() <= 5 );

      return polygonVertices.volume() / geo.volume();
    }


    template< class Geometry, class ReconstructionType, class PointList>
    void computeInterfaceLinePosition ( const Geometry &geo, double concentration, ReconstructionType &g, PointList &intersections )
    {
      intersections.clear();

      double pMin = 0, pMax = 0, volume;
      //use bigger range than [0,1] initially
      double volMin = -1;
      double volMax = 2;

      // Initial guess for p
      for( int i = 0; i < geo.corners(); ++i )
      {
        g.distance() =  -1.0 * ( geo.corner( i ) * g.normal() );

        volume = getVolumeFraction( geo, g );

        if( ( volume <= volMax ) && ( volume >= concentration ) )
        {
          pMax = g.distance();
          volMax = volume;
        }
        if( ( volume >= volMin ) && ( volume <= concentration ) )
        {
          pMin = g.distance();
          volMin = volume;
        }
      }

      g.distance() = brentsMethod( [ &geo, &concentration, &g ] ( double p ) -> double {
                              ReconstructionType h( g.normal(), p );
                              return ( getVolumeFraction( geo, h ) - concentration );
                            }, pMin, pMax, 1e-12 );

      intersections = lineCellIntersections( geo, g );
    }


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRICUTILITY_HH
