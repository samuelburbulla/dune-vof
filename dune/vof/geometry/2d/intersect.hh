#ifndef DUNE_VOF_GEOMETRY_2D_INTERSECT_HH
#define DUNE_VOF_GEOMETRY_2D_INTERSECT_HH

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

#include <dune/vof/geometry/halfspace.hh>
#include <dune/vof/geometry/2d/polygon.hh>

namespace Dune {

  namespace VoF {

    /**
     * \brief namespace containing specific implementations
     */
    namespace __impl {


      /**
       * \ingroup geo2d
       * \brief implementation for an intersection between a polygon and a half space
       *
       * \tparam  Coord  type of the global coordinate
       * \return  the intersection polygon
       */
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


      /**
       * \ingroup geo2d
       * \brief implementation for an intersection between a polygon and a hyperplane
       *
       * \tparam  Coord  type of the global coordinate
       * \return the intersection line
       */
      template< class Coord >
      auto intersect ( const Polygon< Coord >& polygon, const HyperPlane< Coord >& plane ) -> Line< Coord >
      {
        using limits = std::numeric_limits< typename Coord::value_type >;
        using std::abs;

        if ( !plane )
          return Line< Coord >();

        auto container = typename Line< Coord >::Container();

        std::size_t j = 0u;
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

            container[ j++ ] = point ;
          }
          else if ( abs( l0 ) < limits::epsilon() )
            container[ j++ ] = edge.vertex( 0 );
        }

        if ( j < 2 )
          return Line< Coord >();
        else
          return Line< Coord >( std::move( container ) );
      }

    } // namespace __impl

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_2D_INTERSECT_HH
