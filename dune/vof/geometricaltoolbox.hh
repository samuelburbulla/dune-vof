#ifndef __DUNE_GRID_REC_VOL_TRACKEVOLVE_HH__
#define __DUNE_GRID_REC_VOL_TRACKEVOLVE_HH__

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>

namespace Dune
{
  namespace VoF
  {

    template< class V >
    void rotate90degreesCounterClockwise ( V &v )
    {
      std::swap( v[ 0 ], v[ 1 ] );
      v[ 0 ] = -v[ 0 ];
    }

// Gerade der Form 0 = n*x + p
    template< class V >
    class Line2D
    {

    public:

      Line2D< V >()
      {
        n = 0;
        p = 0;
      }

      Line2D< V >( const V &normal, const V &point )
      {
        n = normal;
        p = n * point;
        p = -p;
      }

      Line2D< V >( const Line2D< V >&l )
      {
        n( l.n );
        p( l.p );
      }

      Line2D< V >( const V &n2, const double p2 )
      {
        n = n2;
        p = p2;
      }

      Line2D< V >( const std::pair< V, V > &pair )
      {
        n = pair.first;
        n -= pair.second;
        rotate90degreesCounterClockwise( n );
        p = n * pair.first;
        p = -p;
      }

      V n;
      double p;

    };

    template< class V >
    const V lineIntersection ( const Line2D< V > &g, const Line2D< V > &l  )
    {

      Dune::FieldMatrix< double, 2, 2 > A;
      A[ 0 ][ 0 ] = g.n[ 0 ];
      A[ 0 ][ 1 ] = g.n[ 1 ];
      A[ 1 ][ 0 ] = l.n[ 0 ];
      A[ 1 ][ 1 ] = l.n[ 1 ];

      Dune::FieldVector< double,2  > b;
      b[ 0 ] = -g.p;
      b[ 1 ] = -l.p;


      Dune::FieldVector< double, 2 > x;

      A.solve( x, b );

      return x;
    }


// konvexes 2D-Polygon
    template< class V >
    class Polygon2D
    {

    public:

      Polygon2D< V >() {}

      ~Polygon2D< V >()
      {
        points.clear();
      }

      const V operator[] ( const int i ) const
      {
        return points[ i % points.size() ];
      }


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

            rotate90degreesCounterClockwise( normal );

            if( normal * ( vertex - points[ i ] ) < 0 )
            {
              points.insert( points.begin() + i + 1, vertex );
              return;
            }
          }

      }

      const std::size_t corners () const
      {
        return points.size();
      }

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

        V n = p2 - p1;
        rotate90degreesCounterClockwise( n );

        double skalar = n * ( vertex - p1);

        return skalar >= 0 || std::abs( skalar ) < TOL;

      }

    };


    template< class V >
    void polygonLineIntersection ( const Polygon2D< V > &polygon, const Line2D< V > &g, Polygon2D< V > &intersectionPolygon )
    {
     
      for ( std::size_t i = 0; i < polygon.corners(); ++i )
      {
        if( isOnRecLine( polygon[i], g ) )
        {
          intersectionPolygon.addVertex( polygon[i] );
        }
        else if( isInner( polygon[i], g ) ^ isInner( polygon[i+1], g ) )
        {
          // build line through edge for intersection
          const Line2D< V > lineThroughEdge( std::pair< V, V >( polygon[ i ], polygon[ i+1 ] ) );

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


    template< class GV, class E, class Geo, class V >
    std::vector< V > lineCellIntersections ( const GV &gridView, const E &entity, const Geo &geo, const Line2D< V > &g, const double TOL = 1e-12 )
    {
      const int dim = 2;

      std::vector< V > intersectionPoints;


      const auto &refElement = Dune::ReferenceElements< double, dim >::general( entity.type() );


      for( int k = 0; k < refElement.size( dim-1 ); ++k )
      {
        int i = refElement.subEntity( k, dim-1, 0, dim );
        int j = refElement.subEntity( k, dim-1, 1, dim );

        const V& c0 = geo.global( refElement.position( i, dim ) );
        const V& c1 = geo.global( refElement.position( j, dim ) );


        if( isInner( c0, g, TOL ) ^ isInner( c1, g, TOL ) )
        {

          // build line through edge for intersection
          const Line2D< V > lineThroughEdge( std::pair< V, V >( c0, c1 ) );

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




    template< class V >
    bool isInner ( const V &vertex, const Line2D< V > &g, const double TOL = 1e-12 )
    {
      return vertex * g.n + g.p >= TOL;
    }


    template< class V >
    bool isOnRecLine ( const V &vertex, const Line2D< V > &g, const double TOL = 1e-12 )
    {
      return std::abs( vertex * g.n + g.p ) < TOL;
    }



    template< class Geo, class V >
    void polyAddInnerVertices ( const Geo &geo, const Line2D< V > &g, Polygon2D< V >& polygon, const double TOL = 1e-12 )
    {
      for( int i = 0; i < geo.corners(); ++i )

        if( isInner( geo.corner( i ), g, TOL ) )
          polygon.addVertex( geo.corner( i ) );
    }

    template< class V >
    void polyAddInnerVertices ( const Polygon2D< V >& sourcePolygon, const Line2D< V > &g, Polygon2D< V >& endPolygon, const double TOL = 1e-12 )
    {
      for( std::size_t i = 0; i < sourcePolygon.corners(); ++i )

        if( isInner( sourcePolygon[ i ], g, TOL ) )
          endPolygon.addVertex( sourcePolygon[ i ] );
    }



    template< class GV, class E, class Geo, class V >
    double getVolumeFraction ( const GV &gridView, const E &entity, const Geo &geo, const Line2D< V > &g )
    {
      Polygon2D< V > polygonVertices;


      auto lip = lineCellIntersections( gridView, entity, geo, g );
      for( auto &v : lip )
        polygonVertices.addVertex( v );


      polyAddInnerVertices( geo, g, polygonVertices );


      assert( polygonVertices.corners() <= 5 );

      return polygonVertices.volume() / geo.volume();

    }


  } // end of namespace VoF
} // end of namespace Dune

#endif



