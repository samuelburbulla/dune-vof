#ifndef DUNE_VOF_GEOMETRY_POLYTOPE_HH
#define DUNE_VOF_GEOMETRY_POLYTOPE_HH

#include <functional>
#include <type_traits>

#include <dune/vof/geometry/2d/polygon.hh>
#include <dune/vof/geometry/3d/polyhedron.hh>

namespace Dune
{

  namespace VoF
  {

    // makePolytope
    // ------------
    template < class Geometry >
    static inline auto makePolytope ( const Geometry& geometry )
      -> typename std::enable_if< Geometry::GlobalCoordinate::dimension == 2 && Geometry::mydimension == 2, 
      Polygon< typename Geometry::GlobalCoordinate > >::type
    {
      return makePolygon ( geometry );
    }

    template < class Geometry >
    static inline auto makePolytope ( const Geometry& geometry )
      -> typename std::enable_if< Geometry::GlobalCoordinate::dimension == 2 && Geometry::mydimension == 1, 
      Line< typename Geometry::GlobalCoordinate > >::type
    {
      return makeLine ( geometry );
    }

    template < class Geometry, class Map >
    static inline auto makePolytope ( const Geometry& geometry, Map&& map )
      -> typename std::enable_if< Geometry::GlobalCoordinate::dimension == 2,
      Polygon< typename Geometry::GlobalCoordinate > >::type
    {
      return makePolygon ( geometry, std::forward< Map >( map ) );
    }

    template < class Geometry >
    static inline auto makePolytope ( const Geometry& geometry )
      -> typename std::enable_if< Geometry::GlobalCoordinate::dimension == 3, Polyhedron< typename Geometry::GlobalCoordinate > >::type
    {
      return makePolyhedron( geometry );
    }

    template < class Geometry, class Map >
    static inline auto makePolytope ( const Geometry& geometry, Map&& map )
      -> typename std::enable_if< Geometry::GlobalCoordinate::dimension == 3, Polyhedron< typename Geometry::GlobalCoordinate > >::type
    {
      return makePolyhedron( geometry, std::forward< Map >( map ) );
    }


    // return the interface
    // --------------------
    template< class Entity, class ReconstructionSet >
    auto interface( const Entity &entity, const ReconstructionSet &reconstructions )
    {
      auto polygon = makePolytope( entity.geometry() );
      auto it = intersect( std::cref( polygon ), reconstructions[ entity ].boundary() );
      auto interface = static_cast< typename decltype( it )::Result > ( it );
      return interface;
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_POLYTOPE_HH
