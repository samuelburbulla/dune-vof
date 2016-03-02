#ifndef DUNE_VOF_GEOMETRY_INTERSECT_HH
#define DUNE_VOF_GEOMETRY_INTERSECT_HH

#include <type_traits>
#include <vector>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>

#include <dune/vof/geometry/halfspace.hh>
#include <dune/vof/geometry/polytope.hh>
#include <dune/vof/geometry/utility.hh>


namespace Dune {

  namespace VoF {

    namespace Implementation {

      template< class A, class B >
      struct Intersection
      {
        using Result = void_t< A, B >; // placeholder type

        Intersection ( const A& a, const B& b )
        : a_( a ), b_( b )
        {}

        operator Result ()
        {
          DUNE_THROW( NotImplemented, "Intersetion not yet implemented." );
        }

      private:
        const A& a_;
        const B& b_;
      };

      // polytope / halfspace intersection
      // ...

    } // namespace Implementation

    template< class A, class B >
    auto intersect ( const A& a, const B& b ) -> typename Implementation::Intersection< A, B >
    {
      return Implementation::Intersection< A, B >( a, b );
    }


    // old implemetations
    // ------------------

     // polygon / half space   intersection (incomplete - part 1)
    template< class DomainVector, template <class> class HalfSpace >
    void polygonLineIntersection ( const Polygon2D< DomainVector > &polygon, const HalfSpace< DomainVector > &g, Polygon2D< DomainVector > &intersectionPolygon )
    {
      using std::abs;

      for ( std::size_t i = 0; i < polygon.corners(); ++i )
      {
        auto ci  = g.levelSet( polygon[ i ] );
        auto cip = g.levelSet( polygon[ i+1 ] );

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

        auto c0 = g.levelSet( x0 );
        auto c1 = g.levelSet( x1 );

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
        if( g.levelSet( sourcePolygon[ i ] ) > 0.0 )
          endPolygon.addVertex( sourcePolygon[ i ] );
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_INTERSECT_HH
