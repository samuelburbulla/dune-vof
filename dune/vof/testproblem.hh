#ifndef __DUNE_GRID_REC_VOL_TRACK_TESTPROBLEM_HH__
#define __DUNE_GRID_REC_VOL_TRACK_TESTPROBLEM_HH__

#include <cmath>

namespace Dune
{
  namespace VoF
  {

    //Initial Condition
    template< class DomainVector >
    double f ( const DomainVector &x, const double t )
    {
      DomainVector center { 0.5, 0.5 };
      DomainVector offset { 0, 0.0 };
      DomainVector normal { 0.0, 0 };

      center.axpy( std::cos( ( 2 * M_PI / 10 ) * t ), offset );
      center.axpy( - std::sin( ( 2 * M_PI / 10 ) * t ), normal );

      double dist = ( x - center ).two_norm();

      return (dist < 0.15) ? 1 : 0;	
      //return ( dist < 0.4 && ( x[0] < 0.45 || x[0] > 0.55 || ( x[0] > 0.45 && x[0] < 0.55 && x[1] < 0.5 ) ) ) ? 1 : 0;

    }
   



    //Initial Condition
    template< class V >
    double f0 ( const V &x )
    {
      return f( x, 0.0 );
    }



    template< class DomainVector >
    DomainVector psi ( const DomainVector &x, const double t = 0.0 )
    {

      DomainVector r;

      r[ 0 ] =  x[ 1 ] - 0.5;
      r[ 1 ] =  0.5 - x[ 0 ];

      r *= 2 * M_PI / 10;

      return r;
    }

    static const double psiMax ()
    {
      return 2 * M_PI / 10 * sqrt( 2.0 );
    }


  } // end of namespace Vof
} // end of namespace Dune

#endif







