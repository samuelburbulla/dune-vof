#ifndef DUNE_VOF_GEOMETRY_UTILITY_HH
#define DUNE_VOF_GEOMETRY_UTILITY_HH

#include <cassert>
#include <type_traits>

#include <dune/common/deprecated.hh>
#include <dune/common/typetraits.hh>

namespace Dune
{
  namespace VoF
  {

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


      template< class ... Ts >
      struct isHomogeneous : isHomogeneousHelper< Ts ... >
      {};

      template<>
      struct isHomogeneous<> : std::false_type
      {};


      // GCPImpl
      // -------

      template< class T, std::size_t I, class SFINAE = std::enable_if_t< T::dimension-1 == I > >
      struct GCPImpl;

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

    template< class ... Coords >
    auto generalizedCrossProduct ( const Coords& ... coords ) -> typename __impl::GCPImpl< typename __impl::isHomogeneous< Coords ... >::type, sizeof...(Coords) >::Coordinate
    {
      return __impl::GCPImpl< typename __impl::isHomogeneous< Coords ... >::type, sizeof...(Coords) >::apply( coords ... );
    }


    // normalize
    // ---------

    template< class Coord >
    auto normalize ( Coord &normal ) -> void_t< decltype( std::declval< Coord >().two_norm() ) >
    {
      normal /= normal.two_norm();
    }


  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_UTILITY_HH
