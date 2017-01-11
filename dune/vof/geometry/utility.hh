#ifndef DUNE_VOF_GEOMETRY_UTILITY_HH
#define DUNE_VOF_GEOMETRY_UTILITY_HH

#include <cassert>
#include <type_traits>
#include <utility>

#include <dune/common/deprecated.hh>
#include <dune/common/typetraits.hh>

#include "dune/vof/geometry/2d/polygon.hh"
#include "dune/vof/geometry/3d/polyhedron.hh"

namespace Dune
{
  namespace VoF
  {

    /**
     * /brief namespace containing specific implementations
     */
    namespace __impl {

      // isHomogeneous
      // -------------

      template< class T, class ... Us >
      struct isHomogeneousHelper;

      template< class T >
      struct isHomogeneousHelper< T > : std::true_type
      {
        using type = T;
      };

      template< class T, class ... Us >
      struct isHomogeneousHelper< T, T, Us ... > : isHomogeneousHelper< T, Us ... >
      {};

      template< class T, class U, class ... Vs >
      struct isHomogeneousHelper< T, U, Vs ... > : std::false_type
      {};


      /**
       * \brief   check if all types in a parameter pack are the same
       *
       * \tparam  ...Ts   pack of types
       */
      template< class ... Ts >
      struct isHomogeneous : isHomogeneousHelper< Ts ... >
      {};

      /**
       * \brief   specialisation of isHomogeneous for an empty pack
       */
      template<>
      struct isHomogeneous<> : std::false_type
      {};


      // GCPImpl
      // -------
       /**
        * \brief generalized cross product implementation
        *
        * \tparam   T   type of coordinate
        * \tparam   I   number of arguments
        */
      template< class T, std::size_t I, class SFINAE = std::enable_if_t< T::dimension-1 == I > >
      struct GCPImpl;

      /**
       * \brief generalized cross product implementation for 2D
       */
      template< class T >
      struct GCPImpl< T, 1 >
      {
        static_assert( T::dimension == 2, "" );
        using Coordinate = T;

        static Coordinate apply( const Coordinate& v )
        {
          return { -v[ 1 ], v[ 0 ] };
        }
      };

      /**
       * \brief generalized cross product implementation for 3D
       */
      template< class T >
      struct GCPImpl< T, 2 >
      {
        static_assert( T::dimension == 3, "" );
        using Coordinate = T;

        static Coordinate apply( const Coordinate& v, const Coordinate& w )
        {
          return { v[ 1 ] * w[ 2 ] - v[ 2 ] * w[ 1 ],
                   v[ 2 ] * w[ 0 ] - v[ 0 ] * w[ 2 ],
                   v[ 0 ] * w[ 1 ] - v[ 1 ] * w[ 0 ] };
        }
      };


    } // namespace __impl


    // generalizedCrossProduct
    // -----------------------

    /**
     * \ingroup Geometry
     * \brief generalized cross product
     *
     * \tparam  ...Coords   global coordinates
     * \param   coords      pack of coordinates
     */
    template< class ... Coords >
    auto generalizedCrossProduct ( const Coords& ... coords ) -> typename __impl::GCPImpl< typename __impl::isHomogeneous< Coords ... >::type, sizeof...(Coords) >::Coordinate
    {
      return __impl::GCPImpl< typename __impl::isHomogeneous< Coords ... >::type, sizeof...(Coords) >::apply( coords ... );
    }


    // normalize
    // ---------

    /**
     * \ingroup Geometry
     * \brief normalize vector
     *
     * \tparam  Coord   global coordinate type
     * \param   v       vector
     */
    template< class Coord >
    auto normalize ( Coord &v ) -> void_t< decltype( std::declval< Coord >().two_norm() ) >
    {
      v /= v.two_norm();
    }



    // outerProduct
    // ------------
    template < class ctype, int dim >
    FieldMatrix< ctype, dim, dim > outerProduct ( const FieldVector< ctype, dim > &a, const FieldVector< ctype, dim > &b )
    {
      FieldMatrix< ctype, dim, dim > m( 0.0 );
      for ( std::size_t i = 0; i < dim; ++i )
        m[ i ].axpy( a[ i ], b );

      return m;
    }




    // makePolytope
    // ------------
    template < class Geometry >
    static inline auto makePolytope ( const Geometry& geometry )
      -> typename std::enable_if< Geometry::GlobalCoordinate::dimension == 2, Polygon< typename Geometry::GlobalCoordinate > >::type
    {
      return makePolygon ( geometry );
    }

    template < class Geometry, class Map >
    static inline auto makePolytope ( const Geometry& geometry, Map&& map )
      -> typename std::enable_if< Geometry::GlobalCoordinate::dimension == 2, Polygon< typename Geometry::GlobalCoordinate > >::type
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
     -> typename decltype( intersect( makePolytope( entity.geometry() ), reconstructions[ entity ].boundary() ) )::Result
    {
      auto polygon = makePolytope( entity.geometry() );
      auto it = intersect( polygon, reconstructions[ entity ].boundary() );
      auto interface = static_cast< typename decltype( it )::Result > ( it );
      return interface;
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_UTILITY_HH
