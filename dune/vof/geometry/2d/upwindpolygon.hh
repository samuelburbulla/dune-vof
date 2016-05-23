#ifndef DUNE_VOF_GEOMETRY_2D_UPWINDPOLYGON_HH
#define DUNE_VOF_GEOMETRY_2D_UPWINDPOLYGON_HH

//- dune-common includes
#include <dune/common/fvector.hh>

//- local includes
#include <dune/vof/geometry/2d/polygon.hh>
#include <dune/vof/geometry/utility.hh>

namespace Dune
{
  namespace VoF
  {

      // UpwindPolygon
      // -------------

      /**
       * \brief generate upwind polygon for 2d
       */

    template< class IntersectionGeometry, class Coordinate >
    inline Polygon< Coordinate > upwindPolygon2d ( const IntersectionGeometry& iGeometry, const Coordinate& v )
    {
      if ( ( generalizedCrossProduct( iGeometry.corner( 1 ) - iGeometry.corner( 0 ) ) * v ) < 0.0 )
        return Polygon< Coordinate >( { iGeometry.corner( 0 ), iGeometry.corner( 1 ), iGeometry.corner( 1 ) - v, iGeometry.corner( 0 ) - v } );
      else
        return Polygon< Coordinate >( { iGeometry.corner( 1 ), iGeometry.corner( 0 ), iGeometry.corner( 0 ) - v, iGeometry.corner( 1 ) - v } );
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_2D_UPWINDPOLYGON_HH
