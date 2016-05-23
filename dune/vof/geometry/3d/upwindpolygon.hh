#ifndef DUNE_VOF_GEOMETRY_3D_UPWINDPOLYGON_HH
#define DUNE_VOF_GEOMETRY_3D_UPWINDPOLYGON_HH

//- dune-common includes
#include <dune/common/fvector.hh>

//- local includes
#include <dune/vof/geometry/3d/polyhedron.hh>
#include <dune/vof/geometry/utility.hh>

namespace Dune
{
  namespace VoF
  {

      // UpwindPolygon
      // -------------

      /**
       * \brief generate upwind polygon for 3d
       */

    template< class IntersectionGeometry, class Coordinate >
    inline auto upwindPolygon3d ( const IntersectionGeometry& iGeometry, const Coordinate& v )
    {
      std::vector< Coordinate > nodes;

      if ( generalizedCrossProduct( iGeometry.corner( 1 ) - iGeometry.corner( 0 ), iGeometry.corner( 2 ) - iGeometry.corner( 0 ) ) * v < 0.0 )
      {
        for ( int i = 0; i < iGeometry.corners(); ++i )
        nodes.push_back( iGeometry.corner( i ) - v );

        for ( int i = 0; i < iGeometry.corners(); ++i )
        nodes.push_back( iGeometry.corner( i ) );
      }
      else
      {
        for ( int i = 0; i < iGeometry.corners(); ++i )
          nodes.push_back( iGeometry.corner( i ) );

        for ( int i = 0; i < iGeometry.corners(); ++i )
          nodes.push_back( iGeometry.corner( i ) - v );
      }


      auto type = iGeometry.type();
      if( type.isTriangle() )
      {
        const std::vector< std::array< std::size_t, 2 > > edges {
          {{ 0, 3 }}, {{ 1, 4 }}, {{ 2, 5 }}, {{ 0, 1 }}, {{ 0, 2 }}, {{ 1, 2 }}, {{ 3, 4 }}, {{ 3, 5 }}, {{ 4, 5 }},
          {{ 3, 0 }}, {{ 4, 1 }}, {{ 5, 2 }}, {{ 1, 0 }}, {{ 2, 0 }}, {{ 2, 1 }}, {{ 4, 3 }}, {{ 5, 3 }}, {{ 5, 4 }}
        };

        const std::vector< std::vector< std::size_t > > faces {
          {{ 9, 3, 1, 15 }}, {{ 0, 7, 11, 13 }}, {{ 5, 2, 17, 10 }}, {{ 4, 14, 12 }}, {{ 6, 8, 16 }}
        };
        return Polyhedron< Coordinate >( faces, edges, nodes );
      }
      else if( type.isQuadrilateral() )
      {
        const std::vector< std::array< std::size_t, 2 > > edges {
          {{ 0, 4 }}, {{ 1, 5 }}, {{ 2, 6 }}, {{ 3, 7 }}, {{ 0, 2 }}, {{ 1, 3 }}, {{ 0, 1 }}, {{ 2, 3 }}, {{ 4, 6 }}, {{ 5, 7 }},
          {{ 4, 5 }}, {{ 6, 7 }}, {{ 4, 0 }}, {{ 5, 1 }}, {{ 6, 2 }}, {{ 7, 3 }}, {{ 2, 0 }}, {{ 3, 1 }}, {{ 1, 0 }}, {{ 3, 2 }},
          {{ 6, 4 }}, {{ 7, 5 }}, {{ 5, 4 }}, {{ 7, 6 }}
        };

        const std::vector< std::vector< std::size_t > > faces {
          {{ 0, 8, 14, 16 }}, {{ 5, 3, 21, 13 }}, {{ 6, 1, 22, 12 }}, {{ 19, 2, 11, 15 }}, {{ 4, 7, 17, 18 }}, {{ 10, 9, 23, 20 }}
        };
        return Polyhedron< Coordinate >( faces, edges, nodes );
      }
      else
        DUNE_THROW( InvalidStateException, "Invalid GeometryType." );
    }


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_3D_UPWINDPOLYGON_HH
