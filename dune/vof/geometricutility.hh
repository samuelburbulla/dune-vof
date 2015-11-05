#ifndef DUNE_VOF_GEOMETRICUTILITY_HH
#define DUNE_VOF_GEOMETRICUTILITY_HH

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/fmatrix.hh>
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


    template< template <class> class HyperSurface, class fvector >
    const fvector lineIntersection ( const HyperSurface< fvector > &g, const HyperSurface< fvector > &l )  // make dim-universial
    {

      assert ( fvector::dimension == 2 );

      Dune::FieldMatrix< double, 2, 2 > A( { g.normal(), l.normal() }  );
      fvector b( { -g.p(), -l.p() } );
      fvector x;

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
          // build line through edge for intersection
          V normal = polygon[ i ] - polygon[ i+1 ];
          Dune::VoF::rotate90degreesCounterClockwise<V>( normal );
          const HyperSurface< V > lineThroughEdge( normal, polygon[ i ] );

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


    template< class GV, class E, class Geo, template <class> class HyperSurface, class V >
    std::vector< V > lineCellIntersections ( const GV &gridView, const E &entity, const Geo &geo, const HyperSurface<V> &g, const double TOL = 1e-12 )
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
          V normal = c0 - c1;
          Dune::VoF::rotate90degreesCounterClockwise<V>( normal );
          const HyperSurface< V > lineThroughEdge( normal, c0 );

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
      return vertex * g.normal() + g.p() >= TOL;
    }


    template < template <class> class HyperSurface, class V >
    bool isOnRecLine ( const V &vertex, const HyperSurface<V> &g, const double TOL = 1e-12 )
    {
      return std::abs( vertex * g.normal() + g.p() ) < TOL;
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



    template< class GV, class E, class Geo, template <class> class HyperSurface, class V >
    double getVolumeFraction ( const GV &gridView, const E &entity, const Geo &geo, const HyperSurface<V> &g )
    {
      Polygon2D< V > polygonVertices;


      auto lip = lineCellIntersections( gridView, entity, geo, g );
      for( auto &v : lip )
        polygonVertices.addVertex( v );


      polyAddInnerVertices( geo, g, polygonVertices );


      assert( polygonVertices.corners() <= 5 );

      return polygonVertices.volume() / geo.volume();

    }


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRICUTILITY_HH
