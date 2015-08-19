#ifndef __DUNE_GRID_REC_VOL_TRACK_TESTPROBLEM_HH__
#define __DUNE_GRID_REC_VOL_TRACK_TESTPROBLEM_HH__

#include <cmath>

namespace Dune
{
  namespace VoF
  {
    
	
    template < class ct >
    ct sqr ( const ct x )
    {
      return x * x;
    }

    //Initial Condition
    template < class V >
    double f0 ( const V& x )
    {
      if ( sqr(x[0] - 0.5) + sqr(x[1] - 0.75) <= sqr( 0.15 ) )
	return 1.0;
      else 
	return 0.0;
    }

    //Velocity Field

    /*
    // === Deformation field ===
    template < class ct >
    Dune::FieldVector<double, 2> psi( const Dune::FieldVector <ct, 2>& x, double t )
    {
      Dune::FieldVector<double, 2> r;
      r[0] = sin( 4 * M_PI * (x[0] + 0.5) ) * sin( 4 * M_PI * (x[1] + 0.5) );
      r[1] = cos( 4 * M_PI * (x[0] + 0.5) ) * cos( 4 * M_PI * (x[1] + 0.5) ); 
      
      return r;
    }
    */


    // === Single-Vortex ===
    //Velocity Field
    
    Dune::FieldVector<double,2> psi( const Dune::FieldVector<double,2> &x, double t )
    {
      
      Dune::FieldVector<double,2> r;
      r[0] =  x[1] - 0.5; 	/*- sqr( sin( M_PI * x[0] ) ) * 2 * sin( M_PI * x[1] ) * cos( M_PI * x[1] ); */ 
      r[1] =  0.5 - x[0]; 	/*sqr( sin( M_PI * x[1] ) ) * 2 * sin( M_PI * x[0] ) * cos( M_PI * x[0] ); */
      
      return r;
    }


    double psiMax( )
    {
      return 1.0;
    }


  }
}

#endif







