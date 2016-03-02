#ifndef DUNE_VOF_GEOMETRY_ALGORITHM_HH
#define DUNE_VOF_GEOMETRY_ALGORITHM_HH

#include <cassert>
#include <limits>
#include <vector>

#include <dune/vof/brents.hh>
#include <dune/vof/geometry/halfspace.hh>
#include <dune/vof/geometry/polytope.hh>

namespace Dune {

  namespace VoF {

    template< class Geo, template <class> class HalfSpace, class DomainVector >
    double getVolumeFraction ( const Geo &geo, const HalfSpace< DomainVector > &g )
    {
      Polygon2D< DomainVector > polygon;

      for( int i = 0; i < geo.corners(); ++i )
        if( g.levelSet( geo.corner( i ) ) > 0.0 )
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

#endif // #ifndef DUNE_VOF_GEOMETRY_ALGORITHM_HH
