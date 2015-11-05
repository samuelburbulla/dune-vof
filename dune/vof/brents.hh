#ifndef DUNE_VOF_BRENTS_HH
#define DUNE_VOF_BRENTS_HH

//- C++ includes
#include <cassert>
#include <cstdio>


namespace Dune
{
  namespace VoF
  {


    //Quelle: Wikipedia (Brent's-Verfahren)
    template< class F >
    double brentsMethod ( double a, double b, const F &f )
    {
      const double TOL = 1e-14;

      double fa, fb, fc, c, d, e, p, q, m, s, tol, r;


      fa = f( a );
      fb = f( b );

      if( std::abs( fa ) < TOL )
        return a;
      if( std::abs( fb ) < TOL )
        return b;

      assert( fa * fb < 0 );

      c = a; fc = fa;
      d = b - a; e = d;

      int iter = 0;
      int maxiter = 1000;

      while( iter < maxiter )
      {

        iter++;

        if( fb * fc > 0 )
        {
          c = a;
          fc = fa;
          d = b - a;
          e = d;
        }

        if( std::abs( fc ) < std::abs( fb ) )
        {
          a = b;
          b = c;
          c = a;
          fa = fb;
          fb = fc;
          fc = fa;
        }

        tol = 2.0 * 1e-8 * std::abs( b ) + TOL;
        m = ( c - b ) / 2.0;

        if( std::abs( m ) > tol && std::abs( fb ) > 0 )
        {

          if( std::abs( e ) < tol || std::abs( fa ) <= std::abs( fb ) )
          {
            d = m;
            e = m;
          }
          else
          {
            s = fb / fa;

            if( a == c )
            {
              p = 2.0 * m * s;
              q = 1 - s;
            }
            else
            {
              q = fa / fc;
              r = fb / fc;
              p = s * ( 2 * m * q * ( q - r ) - ( b - a ) * ( r - 1 ) );
              q = ( q - 1 ) * ( r - 1 ) * ( s - 1 );
            }

            if( p > 0 )
              q = -q;
            else
              p = -p;

            s = e;
            e = d;

            if( 2.0 * p < 3.0 * m * q - std::abs( tol * q ) && p < std::abs( s * q / 2.0 ) )
              d = p / q;
            else
            {
              d = m;
              e = m;
            }
          }

          a = b;
          fa = fb;

          if( std::abs( d ) > tol )
            b = b + d;
          else
          {
            if( m > 0 )
              b = b + tol;
            else
              b = b - tol;
          }
        }
        else
          break;

        fb = f( b );
      }

      return b;

    }

   } // namespace VoF

} // namespace Dune

#endif // #ifndef DUNE_VOF_BRENTS_HH
