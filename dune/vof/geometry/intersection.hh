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

        auto container = typename Polygon< Coord >::Container();
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

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_INTERSECT_HH
