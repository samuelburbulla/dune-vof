#ifndef DUNE_VOF_GEOMETRY_INTERSECT_HH
#define DUNE_VOF_GEOMETRY_INTERSECT_HH

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <functional>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>

#include <dune/vof/geometry/halfspace.hh>
#include <dune/vof/geometry/utility.hh>
#include <dune/vof/geometry/2d/intersect.hh>
#include <dune/vof/geometry/3d/intersect.hh>


namespace Dune {

  namespace VoF {

    /**
     * \ingroup Geometry
     * \brief expression template for intersection algorithms
     */
    template< class A, class B >
    class GeometricIntersection
    {
      template< class A_, class B_ >
      static auto apply ( const A_& a, const B_& b ) -> decltype( __impl::intersect( std::declval< A_ >(), std::declval< B_ >() ) )
      {
        return __impl::intersect( a, b );
      }

      template< class A_, class B_ >
      static auto apply ( const A_& a, const B_& b ) -> std::enable_if_t< !std::is_same< A_, B_ >::value, decltype( __impl::intersect( std::declval< B_ >(), std::declval< A_ >() ) ) >
      {
        return __impl::intersect( b, a );
      }

      A a_;
      B b_;

    public:
      GeometricIntersection ( A a, B b )
      : a_( a ), b_( b )
      {}

      using Result = decltype( apply( std::cref( a_ ).get(), std::cref( b_ ).get() ) );

      operator Result ()
      {
        return apply( std::cref( a_ ).get() , std::cref( b_ ).get() );
      }
    };


    /**
     * \ingroup Geometry
     * \brief intersect two geometric bodies
     */
    template< class A, class B >
    auto intersect ( A&& a, B&& b ) -> GeometricIntersection< A, B >
    {
      return GeometricIntersection< A, B >( std::forward< A >( a ), std::forward< B >( b ) );
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_GEOMETRY_INTERSECT_HH
