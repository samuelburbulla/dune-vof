#ifndef DUNE_VOF_GEOMETRY_ALGORITHM_HH
#define DUNE_VOF_GEOMETRY_ALGORITHM_HH

#include <cassert>

#include <functional>
#include <limits>
#include <vector>

#include <dune/vof/brents.hh>
#include <dune/vof/geometry/halfspace.hh>
#include <dune/vof/geometry/intersection.hh>
#include <dune/vof/geometry/polygon.hh>

/*
 * TODO:
 * - implement locateHalfSpace( ... )
 * - ...
 */

namespace Dune {

  namespace VoF {

    template< class Coord >
    double getVolumeFraction ( const Polygon< Coord > &polygon, const HalfSpace< Coord > &halfSpace )
    {
      Polygon< Coord > intersection = intersect( std::cref( polygon ), std::cref( halfSpace ) );

      return intersection.volume() / polygon.volume();
    }


    // locateHalfSpace
    // ---------------


    template< class Coord >
    auto locateHalfSpace ( const Polytope< Coord, Coord::dimension >& polytope, const Coord& normal, double fill ) -> HalfSpace< Coord >
    {
      DUNE_THROW( NotImplemented, "locateHalfSpace( ... ) not yet implemented." );
    }

    template< class Coord >
    auto locateHalfSpace ( const Polygon< Coord >& polygon, const Coord& normal, double fill ) -> HalfSpace< Coord >
    {
      using ctype = typename HalfSpace< Coord >::ctype;
      using limits = std::numeric_limits< ctype >;

      ctype pMin = 0, pMax = 0;
      ctype volume,
            volMin = limits::lowest(), // lower bound
            volMax = limits::max(); // upper bound


      // upper/lower bound for
      for( int i = 0; i < polygon.size(); ++i )
      {
        ctype dist = -( normal * polygon.vertex( i ) );
        HalfSpace< Coord > hs( normal, dist );

        volume = getVolumeFraction( polygon, hs );

        if( ( volume <= volMax ) && ( volume >= fill ) )
        {
          pMax = dist;
          volMax = volume;
        }
        if( ( volume >= volMin ) && ( volume <= fill ) )
        {
          pMin = dist;
          volMin = volume;
        }
      }
      assert( volMin <= volMax );

      if ( volMax == limits::max() )
        return HalfSpace< Coord >( normal, pMin );
      else if ( volMin == limits::lowest() )
        return HalfSpace< Coord >( normal, pMax );
      else if ( volMin == volMax )
        return HalfSpace< Coord >( normal, pMin );
      else
        return HalfSpace< Coord >( normal, brentsMethod( [ &polygon, &normal, &fill ] ( ctype p ) -> ctype {
                                                            HalfSpace< Coord > hs( normal, p );
                                                            return ( getVolumeFraction( polygon, hs ) - fill );
                                                          }, pMin, pMax, 1e-12 ) );
    }


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_ALGORITHM_HH
