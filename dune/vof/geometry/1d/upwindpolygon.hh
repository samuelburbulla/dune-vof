#ifndef DUNE_VOF_GEOMETRY_1D_UPWINDPOLYGON_HH
#define DUNE_VOF_GEOMETRY_1D_UPWINDPOLYGON_HH

//- dune-common includes
#include <dune/common/fvector.hh>

//- local includes
#include <dune/vof/geometry/utility.hh>

namespace Dune
{
  namespace VoF
  {

      // UpwindPolygon
      // -------------

      /**
       * \brief generate upwind polygon for 1d
       */

    template< class IntersectionGeometry, class Coordinate >
    inline Line< Coordinate > upwindPolygon1d ( const IntersectionGeometry& iGeometry, const Coordinate& v )
    {
      return Line< Coordinate >( iGeometry.corner( 0 ), iGeometry.corner( 0 ) - v );
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_2D_UPWINDPOLYGON_HH
