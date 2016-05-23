#ifndef DUNE_VOF_UPWINDPOLYGON_HH
#define DUNE_VOF_UPWINDPOLYGON_HH

//- dune-common includes
#include <dune/common/fvector.hh>

//- local includes
#include <dune/vof/geometry/intersect.hh>
#include <dune/vof/geometry/2d/polygon.hh>
#include <dune/vof/geometry/3d/polyhedron.hh>
#include <dune/vof/geometry/2d/upwindpolygon.hh>
#include <dune/vof/geometry/3d/upwindpolygon.hh>
#include <dune/vof/geometry/utility.hh>


namespace Dune
{
  namespace VoF
  {

    // UpwindPolygon
    // -------------

    /**
     * \brief generate upwind polygon
     *
     * \tparam  IntersectionGeometry
     * \param   iGeometry intersection geometry
     * \param   v         upwind shift
     */
    template< class IntersectionGeometry, class Coordinate >
    inline auto upwindPolygon ( const IntersectionGeometry& iGeometry, const Coordinate& v )
      -> typename std::enable_if< Coordinate::dimension == 2, Polygon< typename IntersectionGeometry::GlobalCoordinate > >::type
    {
      return upwindPolygon2d( iGeometry, v );
    }

    template< class IntersectionGeometry, class Coordinate >
    inline auto upwindPolygon ( const IntersectionGeometry& iGeometry, const Coordinate& v )
      -> typename std::enable_if< Coordinate::dimension == 3, Polyhedron< typename IntersectionGeometry::GlobalCoordinate > >::type
    {
      return upwindPolygon3d( iGeometry, v );
    }


    /**
     * \brief volume of truncated upwind polygon
     */
    template < class P, class Reconstruction >
    inline double truncVolume ( const P& upwind, const Reconstruction& halfSpace )
    {
      P intersection = intersect( std::cref( upwind ), std::cref( halfSpace ) );
      return intersection.volume();
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_UPWINDPOLYGON_HH
