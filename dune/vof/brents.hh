#ifndef DUNE_VOF_BRENTS_HH
#define DUNE_VOF_BRENTS_HH

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <limits>
#include <type_traits>
#include <utility>

namespace Dune
{

  namespace VoF
  {

    // brentsMethod
    // ------------

    /**
     * \ingroup Other
     * \brief Brendt's method
     *
     * \tparam  F             functor type
     * \tparam  X             real type
     * \tparam  Y             real type
     * \param   f             functor (X -> Y)
     * \param   a,b           pair of interval bounds and functor value
     * \param   tolerance     tolerance
     * \param   maxIterations maximum number of iterations
     */
    template< class F, class X, class Y >
    inline static auto brentsMethod ( F f, std::pair< X, Y > a, std::pair< X, Y > b, X tolerance,
                                      std::size_t maxIterations = std::numeric_limits< std::size_t >::max() )
      -> typename std::enable_if< std::is_arithmetic< X >::value && std::is_arithmetic< Y >::value, std::pair< X, Y > >::type
    {
      using std::abs;
      using Z = decltype( a.second * b.second );

      assert( a.second * b.second < Z( 0.0 ) );

      const X half = X( 1 ) / X( 2 );

      std::pair< X, Y > c = a;

      X d = b.first - a.first;
      X e = d;

      for( std::size_t i = 0; i < maxIterations; ++i )
      {
        if( b.second * c.second > Z( 0.0 ) )
        {
          c = a;
          d = (b.first - a.first);
          e = d;
        }

        if( abs( c.second ) < abs( b.second ) )
        {
          a = b;
          b = c;
          c = a;
        }

        X htol = std::numeric_limits< X >::epsilon() * abs( b.first ) + tolerance;
        X tol = X( 2 ) * htol;
        X tm = c.first - b.first;
        X m = half * tm;
        if( (abs( tm ) <= X( 4 )*htol) || (b.second == Y( 0.0 )) )
          return b;

        if( (abs( e ) >= tol) && (abs( a.second ) > abs( b.second )) )
        {
          Y p, q, r;
          Y s = b.second / a.second;

          if( a.first == c.first )
          {
            p = static_cast< Y >( tm ) * s;
            q = Y( 1 ) - s;
          }
          else
          {
            q = a.second / c.second;
            r = b.second / c.second;
            p = s * (static_cast< Y >( tm ) * q * (q - r) - static_cast< X >( b.first - a.first ) * (r - Y( 1 )));
            q = ( q - 1 ) * ( r - 1 ) * ( s - 1 );
          }

          if( p > 0 )
            q = -q;
          else
            p = -p;

          s = e;
          e = d;

          if( (Y( 2 ) * p < Y( 3 ) * m * q - abs( tol * q )) && (p < abs( s * q * half )) )
            d = p / q;
          else
            d = e = m;
        }
        else
          d = e = m;

        a = b;

        if( abs( d ) > tol )
          b.first += d;
        else
          b.first += (m > 0 ? tol : -tol);
        b.second = f( b.first );
      }

      return b;
    }

    /**
     * \ingroup Other
     * \brief Brendt's method
     *
     * \tparam  F             functor type
     * \tparam  T             real type
     * \param   f             functor
     * \param   a,b           interval bounds
     * \param   tolerance     tolerance
     * \param   maxIterations maximum number of iterations
     */
    template< class F, class T >
    inline static auto brentsMethod ( F f, T a, T b, T tolerance,
                                      std::size_t maxIterations = std::numeric_limits< std::size_t >::max() )
      -> typename std::enable_if< std::is_arithmetic< T >::value, T >::type
    {
      using std::abs;

      auto fa = f( a );
      if( abs( fa ) < tolerance )
        return a;

      auto fb = f( b );
      if( abs( fb ) < tolerance )
        return b;

      return brentsMethod( f, std::make_pair( a, fa ), std::make_pair( b, fb ), tolerance, maxIterations ).first;
    }

  } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_BRENTS_HH
