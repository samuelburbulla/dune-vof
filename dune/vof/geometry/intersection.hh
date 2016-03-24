#ifndef DUNE_VOF_GEOMETRY_INTERSECT_HH
#define DUNE_VOF_GEOMETRY_INTERSECT_HH

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <functional>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>

#include <dune/vof/geometry/halfspace.hh>
#include <dune/vof/geometry/polygon.hh>
#include <dune/vof/geometry/polytope.hh>
#include <dune/vof/geometry/utility.hh>

/*
 * TODO:
 * - implement intersection implementations (free functions / callable ?)
 * - ...
 */

namespace Dune {

  namespace VoF {

    namespace __impl {

      // intersect( ... ) implementations

      template< class Coord >
      auto intersect ( const Polytope< Coord, Coord::dimension >& polytope, const HalfSpace< Coord >& halfSpace ) -> Polytope< Coord, Coord::dimension >
      {
        DUNE_THROW( NotImplemented, "__impl::intersect( ... ) not yet implemented." );
      }

      template< class Coord >
      auto intersect ( const Polytope< Coord, Coord::dimension >& polytope, const HyperPlane< Coord >& plane ) -> Polytope< Coord, Coord::dimension-1 >
      {
        DUNE_THROW( NotImplemented, "__impl::intersect( ... ) not yet implemented." );
      }

      // polygon implementations

      template< class Coord >
      auto intersect ( const Polygon< Coord >& polygon, const HalfSpace< Coord >& halfSpace ) -> Polygon< Coord >
      {
        if ( !halfSpace )
          return Polygon< Coord >();

        auto container = typename Line< Coord >::Container();
        container.reserve( polygon.size() );

        for( int i = 0; i < polygon.size( 1 ); ++i )
        {
          auto edge = polygon.edge( i );

          auto l0 = halfSpace.levelSet( edge.vertex( 0 ) );
          auto l1 = halfSpace.levelSet( edge.vertex( 1 ) );

          if( l0 > 0.0 )
            container.push_back( edge.vertex( 0 ) );

          if ( ( l0 > 0.0 ) ^ ( l1 > 0.0 ) )
          {
            Coord point;
            point.axpy( -l1 / ( l0 - l1 ), edge.vertex( 0 ) );
            point.axpy(  l0 / ( l0 - l1 ), edge.vertex( 1 ) );

            container.push_back( point );
          }
        }

        if ( container.empty() )
          return Polygon< Coord >();
        else
          return Polygon< Coord >( std::move( container ) );
      }

      template< class Coord >
      auto intersect ( const Polygon< Coord >& polygon, const HyperPlane< Coord >& plane ) -> Line< Coord >
      {
        using limits = std::numeric_limits< typename Coord::value_type >;
        using std::abs;

        if ( !plane )
          return Line< Coord >();

        auto container = typename Line< Coord >::Container();
        container.reserve( 2u );

        for( int i = 0; i < polygon.size( 1 ); ++i )
        {
          auto edge = polygon.edge( i );

          auto l0 = plane.levelSet( edge.vertex( 0 ) );
          auto l1 = plane.levelSet( edge.vertex( 1 ) );

          if ( ( l0 > 0.0 ) ^ ( l1 > 0.0 ) )
          {
            Coord point;
            point.axpy( -l1 / ( l0 - l1 ), edge.vertex( 0 ) );
            point.axpy(  l0 / ( l0 - l1 ), edge.vertex( 1 ) );

            container.push_back( point );
          }
          else if ( abs( l0 ) < limits::epsilon() )
            container.push_back( edge.vertex( 0 ) );
        }

        if ( container.empty() )
          return Line< Coord >();
        else
          return Line< Coord >( std::move( container ) );
      }

    } // namespace __impl


    template< class A, class B >
    class GeometricIntersection
    {
      template< class A_, class B_ >
      static auto apply ( const A_& a, const B_& b ) -> decltype( __impl::intersect( std::declval< A_ >(), std::declval< B_ >() ) )
      {
        return __impl::intersect( a, b );
      }

      template< class A_, class B_ >
      static auto apply ( const A_& a, const B_& b ) -> std::enable_if_t< !std::is_same< A_, B_ >::value, decltype( __impl::intersect( std::declval< B_ >(), std::declval< A_ >() ) ) >
      {
        return __impl::intersect( b, a );
      }

      A a_;
      B b_;

    public:
      GeometricIntersection ( A a, B b )
      : a_( a ), b_( b )
      {}


      using Result = decltype( apply( std::declval< decltype( std::cref( a_ ).get() ) >(),
                                      std::declval< decltype( std::cref( b_ ).get() ) >() ) );

      operator Result ()
      {
        return apply( std::cref( a_ ).get() , std::cref( b_ ).get() );
      }
    };

    template< class A, class B >
    auto intersect ( A&& a, B&& b ) -> GeometricIntersection< A, B >
    {
      return GeometricIntersection< A, B >( std::forward< A >( a ), std::forward< B >( b ) );
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

     // polygon / half space   intersection (incomplete - part 2)
    template< class DomainVector, class HalfSpace >
    void polyAddInnerVertices ( const Polygon2D< DomainVector > &sourcePolygon, const HalfSpace &g, Polygon2D< DomainVector >& endPolygon, const double TOL = 1e-12 )
    {
      for( std::size_t i = 0; i < sourcePolygon.corners(); ++i )
        if( g.levelSet( sourcePolygon[ i ] ) > 0.0 )
          endPolygon.addVertex( sourcePolygon[ i ] );
    }

     // cell / line intersection
    template< class Geo, template <class> class HalfSpace, class DomainVector >
    void lineCellIntersections ( const Geo &geo,
                                 const HalfSpace< DomainVector > &g,
                                 std::vector< DomainVector > &intersections,
                                 const double TOL = 1e-12 )
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

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_INTERSECT_HH
